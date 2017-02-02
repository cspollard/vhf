{-# OPTIONS_GHC -fplugin GHC.TypeLits.KnownNat.Solver #-}
{-# OPTIONS_GHC -fplugin GHC.TypeLits.Normalise #-}
{-# LANGUAGE DataKinds                 #-}
{-# LANGUAGE FlexibleContexts          #-}
{-# LANGUAGE NoMonomorphismRestriction #-}
{-# LANGUAGE RankNTypes                #-}
{-# LANGUAGE ScopedTypeVariables       #-}
{-# LANGUAGE TypeFamilies              #-}
{-# LANGUAGE TypeOperators             #-}

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
import           Data.Vector  (Vector)
import qualified Data.Vector  as V
import           GHC.TypeLits
import           Linear       as X
import           Linear.V     as X

type M (n :: Nat) (m :: Nat) a = V n (V m a)

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
  reifyDimNat l1 $
    \(Proxy :: Proxy k1) -> reifyDimNat l2 $
      \(Proxy :: Proxy k2) -> do
        (mV :: M k2 k1 a) <- fromVectorM m
        return $ f mV


reifyMatrix1
  :: forall n a r. KnownNat n
  => Vector (Vector a) -> (forall m. KnownNat m => M n m a -> r) -> Maybe r
reifyMatrix1 m f = do
  let l1 = reflectDim (Proxy :: Proxy n)
      l2 = if l1 == 0 then 0 else V.length (m V.! 0)
  reifyDimNat l2 $
    \(Proxy :: Proxy k2) -> do
      (mV :: M n k2 a) <- fromVectorM m
      return $ f mV


reifyMatrix2
  :: forall m a r. KnownNat m
  => Vector (Vector a) -> (forall n. KnownNat n => M n m a -> r) -> Maybe r
reifyMatrix2 m f = do
  let l1 = V.length m
  reifyDimNat l1 $
    \(Proxy :: Proxy k1) -> do
      (mV :: M k1 m a) <- fromVectorM m
      return $ f mV



tailV :: V n a -> V (n-1) a
tailV (V v) = V $ V.tail v

consV :: a -> V n a -> V (n+1) a
consV x (V v) = V $ V.cons x v

headV :: (1 <= n) => V n a -> a
headV (V v) = V.head v

snocV :: V n a -> a -> V (n+1) a
snocV (V v) x = V $ V.snoc v x

catV :: V n a -> V m a -> V (n+m) a
catV (V v) (V v') = V $ v V.++ v'

zeroM :: (Dim m, Dim n, Num a) => M n m a
zeroM = pure zero

{-
-- TODO
-- TODO
-- HERE!
cholesky
  :: forall n a. (Floating a, Conjugate a, KnownNat n, 1 <= n)
  => M n n a -> M n n a
cholesky m =
  where
    r = zero
    s = 1 `consV` r
    t = r `snocV` 1
    -- upp = identity `catV` zeroM
    -- low = transpose $ zeroM `catV` identity
    f :: forall k. (KnownNat k, 1 <= k) => M k k a -> M k k a
    f m =
      if natVal (Proxy :: Proxy n) == 1
        then identity
        else f m
      let t = headV m
          saii = sqrt $ headV t
          bi = (/saii) . conjugate <$> tailV t
          l = saii `consV` bi
          x :: M k (k-1) a
          x = zero `consV` identity
          next :: M (k-1) (k-1) a
          next = transpose . tailV . transpose . tailV $ m
          -- TODO
          -- this is L; need to dot it with the L that comes with (k-1)
      in transpose (consV l (transpose x)) !*! (consV s (transpose (pure 0 `consV` f next)))
-}
