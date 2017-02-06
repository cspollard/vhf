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
import           List.Transformer     (ListT (..), Step (..))
import qualified List.Transformer     as LT
import           Numeric.AD
import           Options.Applicative  hiding (Parser, auto)
import qualified Options.Applicative  as OA
import           System.IO            (IOMode (..), hPutStr, hPutStrLn,
                                       withFile)

import           MarkovChain
import           Metropolis
import           Model
import           Numeric.MCMC


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
    <$> option OA.auto (long "burn")
    <*> option OA.auto (long "skip")
    <*> option OA.auto (long "samples")
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
  c <- next $ LT.drop (n-1) l
  case c of
    Cons x l' -> return . Cons x $ takeEvery n l'
    Nil       -> return Nil


dropWhileL :: Monad m => (a -> Bool) -> ListT m a -> ListT m a
dropWhileL f l = ListT $ do
  c <- next l
  case c of
    Cons x l' ->
      if f x
        then next $ dropWhileL f l'
        else return c
    Nil -> return Nil


parseModel
  :: (FromJSON a, FromJSON b)
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
    Right (dataH :: Vector Int, model :: Model Double, modelparams :: Map Text (ModelParam Double)) -> do

      g <- createSystemRandom

      let s = V.unzip . V.fromList $ M.toList modelparams
          mpnames = fst s
          mps = snd s
          start = fmap _mpInitialValue mps
          logPriors = fmap _mpLogPrior mps
          variations :: Vector (ModelVar Double)
          variations = fmap _mpVariation mps


      let logLH :: forall a. (Floating a, Ord a, Mode a, Scalar a ~ Double) => Vector a -> a
          logLH =
            fromJust
            . modelLogPosterior
                dataH
                (fmap auto model)
                (fmap (fmap auto) variations)
                (fmap (ppToFunc . fmap auto) logPriors)
          gLogLH :: Vector Double -> Vector Double
          gLogLH = grad logLH
          -- xs = take 100 $ conjugateGradientAscent logLH start
          start' = start

      print "start'"
      print start'
      print "llh start"
      print $ logLH start'
      print $ appVars variations start' model


      let c =
            Chain
              (Target logLH $ Just gLogLH)
              (logLH (start' :: Vector Double))
              (start' :: Vector Double)
              Nothing
          nsteps = 10
          eps = 0.01
          trans = hamiltonian eps nsteps
          chain =
            takeEvery nskip
            . LT.drop nburn
            $ runMC trans c g


      withFile outfile WriteMode $ \f -> do
        hPutStrLn f . mconcat . intersperse ", " . fmap T.unpack
          $ "llh" : V.toList mpnames

        LT.runListT . LT.take nsamps
          $ do
            Chain{..} <- chain
            LT.liftIO $ do
              hPutStr f $ show chainScore ++ ", "
              hPutStrLn f
                . mconcat . intersperse ", " . V.toList
                $ show <$> chainPosition
