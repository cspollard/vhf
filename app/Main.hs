{-# LANGUAGE DataKinds                 #-}
{-# LANGUAGE FlexibleContexts          #-}
{-# LANGUAGE NoMonomorphismRestriction #-}
{-# LANGUAGE OverloadedLists           #-}
{-# LANGUAGE OverloadedStrings         #-}
{-# LANGUAGE RecordWildCards           #-}
{-# LANGUAGE ScopedTypeVariables       #-}
{-# LANGUAGE TemplateHaskell           #-}
{-# LANGUAGE TypeFamilies              #-}
{-# LANGUAGE TypeOperators             #-}


module Main where

import           Control.Arrow       ((&&&))
import           Control.Lens
import           Data.List           (intersperse)
import           Data.Map.Strict     (Map)
import qualified Data.Map.Strict     as M
import           Data.Maybe          (fromMaybe)
import           Data.Text           (Text)
import qualified Data.Text           as T
import           GHC.Exts            (IsList (..))
import           GHC.TypeLits
import           Linear
import           Linear.V
import           List.Transformer    (ListT (..), Step (..))
import qualified List.Transformer    as LT
import           Options.Applicative
import           System.IO           (IOMode (..), hPutStr, hPutStrLn, withFile)

import qualified Data.Vector         as V
-- import           Data.Vector.Fixed.Cont (ToPeano (..))

import           MarkovChain
import           Metropolis
import           Probability


type M (n :: Nat) (m :: Nat) a = V m (V n a)
type NC = 4
type NE = 4


fromListV :: forall a n. KnownNat n => [a] -> V n a
fromListV l =
    let v = fromVector $ V.fromList l
        e = error
              $ "failed to convert list to vector of length "
                ++ show (dim (undefined :: V n a))
    in fromMaybe e v

instance KnownNat n => IsList (V n a) where
  type Item (V n a) = a
  fromList = fromListV

zeeSmear :: Floating a => M NC NE a
zeeSmear = transpose . fmap signorm $
  [ [38200, 580, 2.23, 0.0888]
  , [373, 3270, 99.0, 0.851]
  , [4.73, 57.5, 503, 22.5]
  , [0.313, 0.883, 13.6, 101.1]
  ]

zeeSmear2 :: Floating a => M NC NE a
zeeSmear2 = transpose
  [ signorm [0.80, 0.15, 0.0, 0.0]
  , signorm [0.03, 0.90, 0.03, 0.0]
  , signorm [0.0, 0.09, 0.84, 0.07]
  , signorm [0.0, 0.0, 0.02, 0.92]
  ]

zeeData :: V NE Int
zeeData = [44812, 3241, 494, 90]

type Param n m a = (a -> (Model n m a -> Model n m a, a))


-- TODO
-- a model should actually take params as an argument
-- and return a reco spectrum as an output
data Model n m a =
  Model
    { _mBkgs   :: Map Text (V m a)
    , _mSigs   :: V n a
    , _mSmears :: M n m a
    , _mLumi   :: a
    }

makeLenses ''Model


myModel :: Floating a => Model NC NE a
myModel =
  Model
    (M.singleton "ttbar" [1.38e-2, 4.97e-3, 1.20e-3, 5.14e-4])
    (pure 1)
    zeeSmear
    37000


myModelParams :: (Floating a, Ord a) => Map Text (Param NC NE a)
myModelParams = M.fromList
  [ ("ttbarnorm", \x -> (over (mBkgs.ix "ttbar") (fmap (*x)), logLogNormalP 0 0.2 x))
  , ("sigma0", set (mSigs.element 0) &&& nonNegPrior)
  , ("sigma1", set (mSigs.element 1) &&& nonNegPrior)
  , ("sigma2", set (mSigs.element 2) &&& nonNegPrior)
  , ("sigma3", set (mSigs.element 3) &&& nonNegPrior)
  -- , ("smear", \x -> (set mSmears (linearCombM (1-x) x zeeSmear zeeSmear2), logNormalP 0 1 x))
  , ("lumi", \x -> (over mLumi (*x), logLogNormalP 0 0.1 x))
  ]

  where
    nonNegPrior x
      | x < 0 = negate $ 1/0
      | otherwise = 0


myInitialParams :: Fractional a => Map Text a
myInitialParams = M.fromList
  [ ("ttbarnorm", 1)
  , ("sigma0", 1)
  , ("sigma1", 0.5)
  , ("sigma2", 0.1)
  , ("sigma3", 0.01)
  -- , ("smear", 0.0)
  , ("lumi", 1)
  ]

myParamRadii :: Fractional a => Map Text a
myParamRadii = M.fromList
  [ ("ttbarnorm", 0.1)
  , ("sigma0", 0.1)
  , ("sigma1", 0.01)
  , ("sigma2", 0.01)
  , ("sigma3", 0.001)
  -- , ("smear", 0.25)
  , ("lumi", 0.1)
  ]


modelPred :: (KnownNat n, KnownNat m, Num a) => Model n m a -> V m a
modelPred (Model bkgs sigs smears lumi) =
  let bkgTot = foldl (^+^) (pure 0) bkgs
  in fmap (*lumi) $ (+) <$> smears !* sigs <*> bkgTot


appParams :: Num a
          => Map Text (Param n m a)
          -> Model n m a
          -> Map Text a
          -> (Model n m a, a)
appParams fs m ps =
  let fs' = M.intersectionWith ($) fs ps
  in foldr (\(f, p) (model, prior) -> (f model, prior+p)) (m, 0) fs'


modelLogPosterior :: (KnownNat n, KnownNat m, Integral a, Floating b, Ord b)
                  => Map Text (Param n m b)
                  -> Model n m b
                  -> V m a
                  -> Map Text b
                  -> b
modelLogPosterior fs model mdata ps =
  let (model', prior) = appParams fs model ps
      preds = modelPred model'
      logLH = sum $ logPoissonP <$> mdata <*> preds
  in logLH + prior


data InArgs =
    InArgs
        { nburn   :: Int
        , nskip   :: Int
        , nsamps  :: Int
        , outfile :: String
        }

inArgs :: Parser InArgs
inArgs = InArgs
    <$> option auto (long "burn")
    <*> option auto (long "skip")
    <*> option auto (long "samples")
    <*> strOption (long "outfile")

opts :: ParserInfo InArgs
opts = info (helper <*> inArgs) fullDesc

main :: IO ()
main = do

  InArgs {..} <- execParser opts

  {-
  let
      (xs :: [Map Text Double]) = take 100 $ conjugateGradientAscent llh myInitialParams
      start = last xs

  print "data:"
  print zeeData

  print "initial params:"
  mapM_ print $ M.toList myInitialParams

  print "initial prediction:"
  print . modelPred . fst $ appParams myModelParams myModel myInitialParams

  print "gradient ascent:"
  mapM_ print . fmap M.toList $ xs

  print "best fit params:"
  mapM_ print $ M.toList start

  print "best fit prediction:"
  print . modelPred . fst $ appParams myModelParams myModel start
  -}

  let (start :: Map Text Double) = myInitialParams
      llh = modelLogPosterior myModelParams myModel zeeData
      prop = weightedProposal (multibandMetropolis myParamRadii) llh

  g <- createSystemRandom

  let takeEvery n l = ListT $ do
        c <- next $ LT.drop n l
        case c of
          Cons x l' -> return . Cons x $ takeEvery n l'
          Nil       -> return Nil

  let chain = takeEvery nskip . LT.drop nburn $ runMC prop (T start $ llh start) g

  withFile outfile WriteMode $ \f -> do
    hPutStrLn f . mconcat . intersperse ", " . fmap T.unpack $ "llh" : M.keys start

    LT.runListT . LT.take nsamps
      $ do
        (T ps llhxs :: T (Map Text Double) Double) <- chain
        LT.liftIO . hPutStr f $ show llhxs ++ ", "
        LT.liftIO . hPutStrLn f $ mconcat . intersperse ", " . M.elems $ show <$> ps
