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
import           Data.Text            (Text)
import qualified Data.Text            as T
import           Data.Vector          (Vector)
import qualified Data.Vector          as V
import           Linear.Matrix
import           List.Transformer     (ListT (..), Step (..))
import qualified List.Transformer     as LT
import           Numeric.AD
import           Numeric.MCMC
import           Options.Applicative  hiding (Parser, auto)
import qualified Options.Applicative  as OA
import           System.IO            (BufferMode (..), IOMode (..), hPutStr,
                                       hPutStrLn, hSetBuffering, stdout,
                                       withFile)

import           InMatrix             hiding (transpose)
import           MarkovChain
import           Matrix
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
    <$> option OA.auto
      ( long "burn"
      <> help "number of MCMC steps to burn before recording"
      )
    <*> option OA.auto
      ( long "skip"
      <> help "number of MCMC steps to skip between each recording"
      )
    <*> option OA.auto
      ( long "samples"
      <> help "number of samples to record"
      )
    <*> strOption
      ( long "outfile"
      <> help "text file to record to"
      )
    <*> strOption
      ( long "infile"
      <> help "json file to read model from"
      )

opts :: ParserInfo InArgs
opts = info (helper <*> inArgs) fullDesc

main :: IO ()
main = do
  InArgs {..} <- execParser opts

  -- write a line as soon as it comes...
  hSetBuffering stdout LineBuffering

  -- parse the JSON file to Aeson.Values first
  values <- eitherDecode' <$> BS.readFile infile

  -- then try to parsse to our data, Model, and ModelParams
  -- NB: need to give explicit types here so the parser knows what to look for.
  case parseEither parseModel =<< values of
    Left err -> error err
    Right
      ( dataH :: Vector Int
      , model :: Model Double
      , modelparams :: Map Text (ModelParam Double)
      ) -> do

      let (mpnames, mps) = V.unzip . V.fromList $ M.toList modelparams
          start = fmap _mpInitialValue mps
          priors = fmap _mpPrior mps
          variations = fmap _mpVariation mps
          toError = either error id

          -- I'm not sure why we need an explicit type here.
          -- probably because of the RankNType going on here
          logLH
            :: forall a. (Floating a, Ord a, Mode a, Scalar a ~ Double)
            => Vector a -> a
          logLH =
            toError
            . modelLogPosterior
                dataH
                (fmap auto model)
                (fmap (fmap auto) variations)
                (fmap (ppToFunc . fmap auto) priors)

          gLogLH = grad logLH


      -- find the maximum likelihood starting location
      start' <-
        let x = last . take 100 $ conjugateGradientAscent logLH start
        in if any isNaN x || isNaN (logLH x)
          then do
            print "warning: could not find a likelihood maximum"
            print "based on your model and initial parameter values."
            print "using initial values to determine hessian:"
            print "this could be extremely inefficient."
            return start
          else
            return x

      putStrLn ""
      print "starting location:"
      print $ V.zip mpnames start'
      print "starting location likelihood:"
      print $ logLH start'
      putStrLn ""

      -- invert hessian -> covariance matrix
      -- then find the transform from "canonical" variables to "real" variables
      let hess' = hessian (negate . logLH) start'
          cov = toError $ invM hess'
          t = cholM cov
          it = toError $ invM t
          transform v = (t !* v) ^+^ start'
          itransform v' = it !* (v' ^-^ start')

      -- need an RNG...
      g <- createSystemRandom

      -- finally, build the chain, metropolis transition, and the MCMC walk
      let c =
            Chain
              (Target (logLH . transform) $ Just (gLogLH . transform))
              (logLH start')
              (itransform start')
              Nothing

          trans = metropolis (1 / fromIntegral nskip)
          walk = takeEvery nskip . LT.drop nburn $ runMC trans c g


      -- write the walk locations to file.
      withFile outfile WriteMode $ \f -> do
        hPutStrLn f . mconcat . intersperse ", " . fmap T.unpack
          $ "llh" : V.toList mpnames

        LT.runListT . LT.take nsamps $ do
          Chain{..} <- walk
          LT.liftIO $ do
            hPutStr f $ show chainScore ++ ", "
            hPutStrLn f
              . mconcat . intersperse ", " . V.toList
              $ show <$> transform chainPosition


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
