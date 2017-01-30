{-# LANGUAGE DataKinds                 #-}
{-# LANGUAGE NoMonomorphismRestriction #-}
{-# LANGUAGE OverloadedLists           #-}
{-# LANGUAGE OverloadedStrings         #-}
{-# LANGUAGE RecordWildCards           #-}
{-# LANGUAGE ScopedTypeVariables       #-}
{-# LANGUAGE TemplateHaskell           #-}
{-# LANGUAGE TypeFamilies              #-}


module Main where

import           Control.Applicative (liftA2)
import           Control.Arrow       ((&&&))
import           Control.Lens
import           Data.List           (intersperse)
import           Data.Map.Strict     (Map)
import qualified Data.Map.Strict     as M
import           Data.Text           (Text)
import qualified Data.Text           as T
import           Data.Traversable    (mapAccumL)
import           List.Transformer    (ListT (..), Step (..))
import qualified List.Transformer    as LT
import           Options.Applicative
import           System.IO           (IOMode (..), hPutStr, hPutStrLn, withFile)

import           MarkovChain
import           Matrix
import           Metropolis
import           Probability



-- test case

type NE = ToPeano 4
type NC = ToPeano 4


cutOffNormal :: (InvErf b, Variate b, PrimMonad m, Ord b) => b -> b -> Prob m b
cutOffNormal mu s = do
  x <- normal mu s
  if x < 0 then cutOffNormal mu s else return x

zeeSmear :: Fractional a => Mat NC NE a
zeeSmear = transpose
  [ normalize [38200, 580, 2.23, 0.0888]
  , normalize [373, 3270, 99.0, 0.851]
  , normalize [4.73, 57.5, 503, 22.5]
  , normalize [0.313, 0.883, 13.6, 101.1]
  ]

zeeSmear2 :: Fractional a => Mat NC NE a
zeeSmear2 = transpose
  [ normalize [0.80, 0.15, 0.0, 0.0]
  , normalize [0.03, 0.90, 0.03, 0.0]
  , normalize [0.0, 0.09, 0.84, 0.07]
  , normalize [0.0, 0.0, 0.02, 0.92]
  ]

normalize :: (Traversable t, Fractional c) => t c -> t c
normalize v = let (s', v') = mapAccumL (\s x -> (s+x, x/s')) 0 v in v'

zeeData :: Vec NE Int
zeeData = [44812, 3241, 494, 90]

nE :: Int
nE = arity (undefined :: NE)

nC :: Int
nC = arity (undefined :: NC)

type Hist = Vec
type Param n m a = (a -> (Model n m a -> Model n m a, a))


-- TODO
-- a model should actually take params as an argument
-- and return a reco spectrum as an output
data Model n m a =
  Model
    { _mBkgs   :: Map Text (Hist m a)
    , _mSigs   :: Hist n a
    , _mSmears :: Mat n m a
    , _mLumi   :: a
    }

makeLenses ''Model


myModel :: Fractional a => Model NC NE a
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


modelPred :: (Arity n, Arity m, Num a) => Model n m a -> Hist m a
modelPred (Model bkgs sigs smears lumi) =
  let bkgTot = foldl (liftA2 (+)) (pure 0) bkgs
  in fmap (*lumi) $ (+) <$> multMV smears sigs <*> bkgTot


appParams :: Num a
          => Map Text (Param n m a)
          -> Model n m a
          -> Map Text a
          -> (Model n m a, a)
appParams fs m ps =
  let fs' = M.intersectionWith ($) fs ps
  in foldr (\(f, p) (model, prior) -> (f model, prior+p)) (m, 0) fs'


modelLogPosterior :: (Arity m, Arity n, Integral a, Floating b, Ord b)
                  => Map Text (Param n m b)
                  -> Model n m b
                  -> Hist m a
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
