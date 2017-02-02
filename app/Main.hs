{-# LANGUAGE DataKinds                 #-}
{-# LANGUAGE FlexibleContexts          #-}
{-# LANGUAGE NoMonomorphismRestriction #-}
{-# LANGUAGE OverloadedLists           #-}
{-# LANGUAGE OverloadedStrings         #-}
{-# LANGUAGE RankNTypes                #-}
{-# LANGUAGE RecordWildCards           #-}
{-# LANGUAGE ScopedTypeVariables       #-}
{-# LANGUAGE TypeFamilies              #-}
{-# LANGUAGE TypeOperators             #-}


module Main where

import           Data.Aeson
import           Data.Aeson.Types     (Parser, parseEither)
import qualified Data.ByteString.Lazy as BS
import           Data.List            (intersperse)
import           Data.Map.Strict      (Map)
import qualified Data.Map.Strict      as M
import           Data.Maybe           (fromJust)
import           Data.Text            (Text)
import qualified Data.Text            as T
import           Data.Vector          (Vector)
import qualified Data.Vector          as V
import           GHC.Exts             (IsList (..))
import           List.Transformer     (ListT (..), Step (..))
import qualified List.Transformer     as LT
import           Numeric.MCMC
import           Options.Applicative  hiding (Parser)
import qualified Options.Applicative  as OA
import           System.IO            (IOMode (..), hPutStr, hPutStrLn,
                                       withFile)

import           MarkovChain
import           Model


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

everyLT :: Monad m => Int -> (a -> m ()) -> ListT m a -> ListT m a
everyLT n f = go 0
  where
    go m ll = ListT $ do
      c <- next ll
      case c of
        Cons x lll ->
          if m >= n
            then f x >> (return . Cons x $ go 0 lll)
            else return . Cons x $ go (m+1) lll
        Nil -> return Nil


takeEvery :: Monad m => Int -> ListT m a -> ListT m a
takeEvery n l = ListT $ do
  c <- next $ LT.drop n l
  case c of
    Cons x l' -> return . Cons x $ takeEvery n l'
    Nil       -> return Nil


parseModel
  :: (Floating a, FromJSON a, FromJSON b)
  => Value -> Parser (Vector b, Model a, Map Text (ModelParam a))
parseModel = withObject "error: parseModel was not given a json object" $
  \o -> do
    d <- o .: "Data"
    m <- o .: "Nominal"
    mps <- o .: "ModelVars"

    return (d, m, mps)

main :: IO ()
main = do
  InArgs {..} <- execParser opts

  values <- eitherDecode' <$> BS.readFile infile

  case parseEither parseModel =<< values of
    Left err -> print err
    Right (dataH :: Vector Int, model :: Model Double, modelparams) -> do

      g <- createSystemRandom

      let (mpnames, mps) = V.unzip . V.fromList $ M.toList modelparams
          start = mpInitialValue <$> mps
          logPriors = unPP . mpLogPrior <$> mps
          variations = mpVariation <$> mps
          logLH = modelLogPosterior dataH model variations logPriors
          -- sd = 2.84*2.84 / fromIntegral (length start)
          trans = metropolis 0.001
            -- adaptiveMetropolis sd 0.0001 lLH
            -- trans = concatAllT $ replicate nburn (metropolis 0.1) ++ repeat (slice 0.02)
            -- trans = concatAllT $ replicate nburn (metropolis 0.001)
            -- trans = hamiltonian 0.01 5
          c =
            Chain
              (Target (fromJust . logLH) Nothing)
              (fromJust . logLH $ start)
              start
              Nothing
          -- c = T (logLH start) (AMInfo 2 start start start $ outer start start)
          chain =
            takeEvery nskip . LT.drop nburn $ runMC trans c g

      withFile outfile WriteMode $ \f -> do
        hPutStrLn f . mconcat . intersperse ", " . fmap T.unpack
          $ "llh" : V.toList mpnames

        LT.runListT . LT.take nsamps
          $ do
            Chain {..} <- chain
            LT.liftIO . hPutStr f $ show chainScore ++ ", "
            LT.liftIO . hPutStrLn f
              . mconcat . intersperse ", " . toList
              $ show <$> chainPosition
