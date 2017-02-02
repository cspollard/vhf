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
import           Linear
import           List.Transformer     (ListT (..), Step (..))
import qualified List.Transformer     as LT
-- import           Numeric.MCMC
import           Options.Applicative  hiding (Parser)
import qualified Options.Applicative  as OA
import           System.IO            (IOMode (..), hPutStr, hPutStrLn,
                                       withFile)

import           MarkovChain
import           Metropolis
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
          logLH = fromJust . modelLogPosterior dataH model variations logPriors

      let radii = const 0.5 <$> start
          cov0 = outer radii radii
          ami = AMInfo 1 start start start cov0
          c = T (logLH start) ami
          sd = 2.4*2.4 / fromIntegral (length start)
          eps = 0.1
          trans = adaptiveMetropolis (fromIntegral nburn) cov0 sd eps logLH
          chain =
            takeEvery nskip
            . dropWhileL ((< fromIntegral nburn) . amt . sndT)
            $ runMC trans c g


      withFile outfile WriteMode $ \f -> do
        hPutStrLn f . mconcat . intersperse ", " . fmap T.unpack
          $ "llh" : V.toList mpnames

        LT.runListT . LT.take nsamps
          $ do
            T x AMInfo{..} <- chain
            LT.liftIO $ do
              print "t:"
              print amt

              print "pos:"
              print ampost

              print "avg:"
              print amavgt

              print "cov:"
              print amcovt

              hPutStr f $ show x ++ ", "
              hPutStrLn f
                . mconcat . intersperse ", " . V.toList
                $ show <$> ampost
