{-# LANGUAGE DataKinds                 #-}
{-# LANGUAGE DeriveGeneric             #-}
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

import           Control.Arrow           ((&&&))
import           Control.Lens
import           Data.Aeson
import           Data.Aeson.TH
import           Data.Aeson.Types        (Parser, parseEither, typeMismatch)
import qualified Data.ByteString.Lazy    as BS
import           Data.List               (intersperse)
import           Data.Map.Strict         (Map)
import qualified Data.Map.Strict         as M
import           Data.Maybe              (fromMaybe, isNothing)
import           Data.Text               (Text)
import qualified Data.Text               as T
import           GHC.Exts                (IsList (..))
import           GHC.Generics
import           Linear
import           List.Transformer        (ListT (..), Step (..))
import qualified List.Transformer        as LT
import           Numeric.MCMC
import           Options.Applicative     hiding (Parser)
import qualified Options.Applicative     as OA
import           System.IO               (IOMode (..), hPutStr, hPutStrLn,
                                          withFile)

import           MarkovChain
-- import           Metropolis
import           Hamiltonian
import           Numeric.AD              (grad)
import qualified Numeric.AD.Mode.Reverse as R
import           Probability


type Hist = ZipList
type Mat a = Hist (Hist a)

instance IsList (ZipList a) where
  type Item (ZipList a) = a
  fromList = ZipList
  toList = getZipList

instance FromJSON a => FromJSON (ZipList a) where
  parseJSON v = ZipList <$> parseJSON v

instance ToJSON a => ToJSON (ZipList a) where
  toJSON = toJSON . getZipList

type Param a = (a -> (Model a -> Model a, a))


-- TODO
-- a model should actually take params as an argument
-- and return a reco spectrum as an output?
data Model a =
  Model
    { _mBackgrounds :: Map Text (Hist a)
    , _mSignal      :: Hist a
    , _mMigration   :: Mat a
    , _mLuminosity  :: a
    } deriving (Generic, Show)

makeLenses ''Model

instance FromJSON a => FromJSON (Model a) where
  parseJSON = genericParseJSON $ defaultOptions{fieldLabelModifier=drop 2}

instance ToJSON a => ToJSON (Model a) where
  toJSON = genericToJSON $ defaultOptions{fieldLabelModifier=drop 2}


newtype ParamPrior a = ParamPrior { _unPP :: a -> a }

instance (Floating a, FromJSON a) => FromJSON (ParamPrior a) where
  parseJSON (String "Flat") = return . ParamPrior $ const 0

  parseJSON (Object o) =
    (o .: "Normal" >>= p logNormalP)
    <|> (o .: "LogNormal" >>= p logLogNormalP)
    where
      p f (Object m) = fmap ParamPrior $ f <$> m .: "Mu" <*> m .: "Sigma"
      p _ invalid    = typeMismatch "ParamPrior" invalid

  parseJSON invalid = typeMismatch "ParamPrior" invalid


newtype ModelVariation a = ModelVariation { _unMV :: a -> Model a -> Model a }


parseModeling :: (FromJSON a, Num a) => Object -> Parser (ModelVariation a)
parseModeling obj = do
  mmVar <- obj .:? "MigrationMatrix"
  bkgVars <- obj .:? "Backgrounds"


  let
    f x (Model bkgs sig mm lumi) =
      let bkgs' = bkgVars <&> M.unionWith (\v v' -> (((1-x) *^ v) ^+^ (x *^ v'))) bkgs
          mm' = mmVar <&> \m -> ((1-x) *!! mm) !+! (x *!! m)
      in Model
          (fromMaybe bkgs bkgs')
          sig
          (fromMaybe mm mm')
          lumi

  if isNothing mmVar && isNothing bkgVars
    then fail "MigrationMatrix and Backgrounds are missing"
    else return $ ModelVariation f

parseSignal :: (FromJSON a, Num a) => Object -> Parser (ModelVariation a)
parseSignal obj = do
  s <- obj .: "Signal"
  return . ModelVariation $ \x -> over mSignal ((x *^ s) ^+^)


instance (Floating a, FromJSON a) => FromJSON (ModelVariation a) where
  parseJSON (String "Lumi") = return . ModelVariation $ \x -> over mLuminosity (*x)

  parseJSON (Object obj) =
    parseModeling obj <|> parseSignal obj

  parseJSON invalid = typeMismatch "ModelVariation" invalid


data ModelParam a =
  ModelParam
    { _mpInitialValue :: a
    , _mpRadius       :: a
    , _mpPrior        :: ParamPrior a
    , _mpVariation    :: ModelVariation a
    } deriving Generic

makeLenses ''ModelParam

instance (Floating a, FromJSON a) => FromJSON (ModelParam a) where
  parseJSON = genericParseJSON $ defaultOptions{fieldLabelModifier=drop 3}


parseModel :: (Floating a, FromJSON a)
           => Value -> Parser (Hist Int, Model a, Map Text (ModelParam a))
parseModel = withObject "error: decodeModel was not given a json dictionary" $
  \o -> do
    m <- o .: "Model"
    mps <- o .: "ModelParams"
    d <- o .: "Data"

    return (d, m, mps)


-- TODO
-- with ziplists (^+^) is not the same as liftA2 (+)
modelPred :: (Num a) => Model a -> Hist a
modelPred (Model bkgs sigs smears lumi) =
  let bkgTot = foldl (liftA2 (+)) (pure 0) bkgs
  in fmap (*lumi) $ (smears !* sigs) ^+^ bkgTot


appParams :: Num a
          => Map Text (Param a)
          -> Model a
          -> Map Text a
          -> (Model a, a)
appParams fs m ps =
  let fs' = M.intersectionWith ($) fs ps
  in foldr (\(f, p) (model, prior) -> (f model, prior+p)) (m, 0) fs'


modelLogPosterior :: (Integral a, Floating b, Ord b)
                  => Map Text (Param b)
                  -> Model b
                  -> Hist a
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
        , infile  :: String
        }

inArgs :: OA.Parser InArgs
inArgs = InArgs
    <$> option auto (long "burn")
    <*> option auto (long "skip")
    <*> option auto (long "samples")
    <*> strOption (long "outfile")
    <*> strOption (long "infile")

opts :: ParserInfo InArgs
opts = info (helper <*> inArgs) fullDesc

main :: IO ()
main = do

  InArgs {..} <- execParser opts
  parsed <- eitherDecode' <$> BS.readFile infile
  case parseEither parseModel =<< parsed of
    Left err -> print err
    Right (dataH, model, modelparams) -> do

      let start = _mpInitialValue <$> modelparams
          radii = _mpRadius <$> modelparams
          params =
            (\p -> (_unMV . _mpVariation) p &&& (_unPP . _mpPrior) p) <$> modelparams
          llh = modelLogPosterior params model dataH

      let trans = slice 0.1

      g <- createSystemRandom

      let takeEvery n l = ListT $ do
            c <- next $ LT.drop n l
            case c of
              Cons x l' -> return . Cons x $ takeEvery n l'
              Nil       -> return Nil

      let c = Chain (Target llh Nothing) (llh start) start Nothing
          chain = takeEvery nskip . LT.drop nburn $ runMC trans c g

      withFile outfile WriteMode $ \f -> do
        hPutStrLn f . mconcat . intersperse ", " . fmap T.unpack $ "llh" : M.keys start

        LT.runListT . LT.take nsamps
          $ do
            Chain{..} <- chain
            LT.liftIO . print . modelPred . fst . appParams params model $ chainPosition
            LT.liftIO . hPutStr f $ show chainScore ++ ", "
            LT.liftIO . hPutStrLn f $ mconcat . intersperse ", " . M.elems $ show <$> chainPosition



{-
myModel :: Floating a => Model a
myModel =
  Model
    (M.singleton "ttbar" [1.38e-2, 4.97e-3, 1.20e-3, 5.14e-4])
    (pure 1)
    zeeSmear
    37000


myModelParams :: (Floating a, Ord a) => Map Text (Param a)
myModelParams = M.fromList
  [ ("ttbarnorm", \x -> (over (mBackgrounds.ix "ttbar") (fmap (*x)), logLogNormalP 0 0.2 x))
  , ("sigma0", set (mSignal.element 0) &&& nonNegPrior)
  , ("sigma1", set (mSignal.element 1) &&& nonNegPrior)
  , ("sigma2", set (mSignal.element 2) &&& nonNegPrior)
  , ("sigma3", set (mSignal.element 3) &&& nonNegPrior)
  -- , ("smear", \x -> (set mSmears (linearCombM (1-x) x zeeSmear zeeSmear2), logNormalP 0 1 x))
  , ("lumi", \x -> (over mLuminosity (*x), logLogNormalP 0 0.1 x))
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

zeeSmear :: Floating a => Mat a
zeeSmear = sequenceA . fmap signorm $
  [ [38200, 580, 2.23, 0.0888]
  , [373, 3270, 99.0, 0.851]
  , [4.73, 57.5, 503, 22.5]
  , [0.313, 0.883, 13.6, 101.1]
  ]

zeeSmear2 :: Floating a => Mat a
zeeSmear2 = sequenceA . fmap signorm $
  [ [0.80, 0.15, 0.0, 0.0]
  , [0.03, 0.90, 0.03, 0.0]
  , [0.0, 0.09, 0.84, 0.07]
  , [0.0, 0.0, 0.02, 0.92]
  ]

zeeData :: Hist Int
zeeData = [44812, 3241, 494, 90]

-}
