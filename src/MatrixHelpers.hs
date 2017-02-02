{-# LANGUAGE NoMonomorphismRestriction #-}
{-# LANGUAGE RankNTypes                #-}
{-# LANGUAGE ScopedTypeVariables       #-}

module MatrixHelpers
  ( module X
  , M
  , toVectorM
  , fromVectorM
  , reifyMatrix
  , reifyMatrix1
  , reifyMatrix2
  ) where

import           Data.Proxy
import           Data.Vector (Vector)
import qualified Data.Vector as V
import           Linear      as X
import           Linear.V    as X

type M n m a = V m (V n a)

fromVectorM :: (Dim n, Dim m) => Vector (Vector a) -> Maybe (M n m a)
fromVectorM m = do
  vs <- traverse fromVector m
  fromVector vs

toVectorM :: M n m a -> Vector (Vector a)
toVectorM = toVector . fmap toVector


reifyMatrix :: Vector (Vector a) -> (forall n m. (Dim n, Dim m) => M n m a -> r) -> Maybe r
reifyMatrix m f = do
  let l1 = V.length m
      l2 = if l1 == 0 then 0 else V.length (m V.! 0)
  reifyDim l1 $
    \(Proxy :: Proxy k1) -> reifyDim l2 $
      \(Proxy :: Proxy k2) -> do
        (mV :: M k2 k1 a) <- fromVectorM m
        return $ f mV


reifyMatrix1
  :: forall n a r. Dim n
  => Vector (Vector a) -> (forall m. Dim m => M n m a -> r) -> Maybe r
reifyMatrix1 m f = do
  let l1 = reflectDim (Proxy :: Proxy n)
      l2 = if l1 == 0 then 0 else V.length (m V.! 0)
  reifyDim l2 $
    \(Proxy :: Proxy k2) -> do
      (mV :: M n k2 a) <- fromVectorM m
      return $ f mV


reifyMatrix2
  :: forall m a r. Dim m
  => Vector (Vector a) -> (forall n. Dim m => M n m a -> r) -> Maybe r
reifyMatrix2 m f = do
  let l1 = V.length m
  reifyDim l1 $
    \(Proxy :: Proxy k1) -> do
      (mV :: M k1 m a) <- fromVectorM m
      return $ f mV
