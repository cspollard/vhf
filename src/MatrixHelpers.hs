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
  , cholesky
  ) where

import           Control.Applicative
import           Control.Lens
import           Data.Monoid
import           Data.Proxy
import           Data.Vector         (Vector)
import qualified Data.Vector         as V
import           GHC.TypeLits
import           Linear              as X
import           Linear.Trace        as X
import           Linear.V            as X

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



cholesky
  :: forall n a. (Floating a, Conjugate a, KnownNat n, 1 <= n)
  => M n n a -> Maybe (M n n a)
cholesky a = l
  where
    l = f a

    l2 = liftA2 (!*!) l l
    midx m i j = (m ^? ix i) >>= (^? ix j)

    -- TODO
    -- why is the trace instance for V failing?
    -- da = diagonal a
    -- dl = diagonal <$> l

    f :: M n n a -> Maybe (M n n a)
    f = itraverse (itraverse . h)

    h :: Int -> Int -> a -> Maybe a
    h i j _
      | i < j = Just 0

      | i == j = do
        x <- midx a i i
        y <- s i j
        return . sqrt $ x - y

      | i > j = do
        l' <- l
        x <- midx l' i i
        y <- midx a i j
        z <- s i j
        return $ (y - z) / x


    s i j = do
      v <- (^? ix i) =<< l2
      return . getSum $ foldMapOf (itakingWhile (\k _ -> k < j) folded) Sum v


{-
-- TODO
-- TODO
-- this is not at all validated.
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
      in transpose (consV l (transpose x)) !*! (consV s (transpose (pure 0 `consV` f next)))
-}
