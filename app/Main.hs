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

import           Control.Arrow        ((&&&))
import           Control.Lens
import           Data.Aeson
import           Data.Aeson.TH
import           Data.Aeson.Types     (Parser, parseEither, typeMismatch)
import qualified Data.ByteString.Lazy as BS
import qualified Data.Foldable        as F
import           Data.List            (intersperse)
import           Data.Map.Strict      (Map)
import qualified Data.Map.Strict      as M
import           Data.Maybe           (fromMaybe, isNothing)
import           Data.Text            (Text)
import qualified Data.Text            as T
import qualified Data.Vector          as V
import           GHC.Exts             (IsList (..))
import           GHC.Generics
import           GHC.TypeLits
import           Linear
import           Linear.V
import           List.Transformer     (ListT (..), Step (..))
import qualified List.Transformer     as LT
import           Numeric.AD           hiding (auto)
import           Numeric.MCMC
import           Options.Applicative  hiding (Parser)
import qualified Options.Applicative  as OA
import           System.IO            (BufferMode (..), IOMode (..), hPutStr,
                                       hPutStrLn, hSetBuffering, stdout,
                                       withFile)

import           MarkovChain
import           Metropolis
import           Model
import           Probability


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

main :: IO ()
main = do
  InArgs {..} <- execParser opts

  parsed <- eitherDecode' <$> BS.readFile infile

  case parseEither parseModel =<< parsed of
    Left err -> print err
    Right (dataH, model, modelparams) -> do

      g <- createSystemRandom

      let keys = M.keys modelparams
          legend = zip keys [0..]
          mps = fromList $ M.elems modelparams
          start = _mpInitialValue <$> mps :: V n Double
          params = mps <&>
            (\p -> (_unMV . _mpVariation) p &&& (_unPP . _mpPrior) p)
          llh = modelLogPosterior params model dataH
          sd = 2.84*2.84 / fromIntegral (length start)
          trans = adaptiveMetropolis sd 0.0001 llh
            -- metropolis 0.001
            -- trans = concatAllT $ replicate nburn (metropolis 0.1) ++ repeat (slice 0.02)
            -- trans = concatAllT $ replicate nburn (metropolis 0.001)
            -- trans = hamiltonian 0.01 5
          -- c = Chain (Target llh Nothing) (llh start) start Nothing
          c = T (llh start) (AMInfo 2 start start start $ outer start start)
          chain =
            takeEvery nskip . LT.drop nburn $ runMC trans c g

      withFile outfile WriteMode $ \f -> do
        hPutStrLn f . mconcat . intersperse ", " . fmap T.unpack $ "llh" : keys

        LT.runListT . LT.take nsamps
          $ do
            T x AMInfo{..} <- chain
            -- LT.liftIO . print . modelPred . fst . appParams params model $ chainPosition
            LT.liftIO . hPutStr f $ show x ++ ", "
            LT.liftIO . hPutStrLn f
              . mconcat . intersperse ", " . toList
              $ show <$> ampost
