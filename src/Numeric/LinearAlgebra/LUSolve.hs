--
-- LUSolve.hs
--
-- Solve an LU factored system of equations.
--

module Numeric.LinearAlgebra.LUSolve where

import           Control.Loop                      (numLoop, numLoopState)
import           Control.Monad.ST                  (ST, runST)
import qualified Data.Matrix.Dense.Generic         as M
import qualified Data.Matrix.Dense.Generic.Mutable as MU
import qualified Data.Vector.Unboxed               as V
import qualified Data.Vector.Unboxed.Mutable       as VU



-- Factor a rectangular matrix A into PA = LU, where L is lower triangular,
-- U is upper triangular and P is the permuation matrix (represented as
-- a vector containing the location of the nonzero column) that describes
-- the row interchanges (partial pivoting).
--
luFactor :: M.Matrix V.Vector Double -> (M.Matrix V.Vector Double, V.Vector Int)
luFactor aOrig = runST $ do
    -- the first thing to do is to copy the original immutable
    -- matrix into a mutable one, since the LU factorization runs
    -- efficiently in place. At the end of the calculation, the
    -- matrix will be frozen (i.e., marked as immutable again).

    let
        (nr, nc) = M.dim aOrig

    a <- M.thaw aOrig
    p <- V.unsafeThaw $ V.generate nr id  -- initialize the permutation vector
                                          -- to the identity, [0,1,2,...,(nr - 1)]

    --luFactor_ a 0 (nr - 1) 0 (nc - 1) p

    a' <- M.unsafeFreeze a
    p' <- V.unsafeFreeze p

    return (a', p')


luFactor_ :: MU.MMatrix VU.MVector s Double
          -> Int
          -> Int
          -> Int
          -> Int
          -> VU.MVector s Int -> ST s ()
luFactor_ a rl rh cl ch p = do
    luFactor_ a rl rh cl ch p
    rowSwap a 1 2 undefined Up
    triangularSolve a undefined undefined
    matrixMultiply 1.0 (subMatrix (0,0) (1,1) a)
                       (subMatrix (0,0) (1,1) a)
                   1.0 (subMatrix (0,0) (1,1) a)
    luFactor_ a rl rh cl ch p
    rowSwap a 1 2 undefined Up


-- This is a generic mutable matrix multiply.  The mutable references
-- need not be distinct, allowing it to be used as part of an in-place
-- algorithm like LU factorization.
--
-- What is computed (and stored in the matrix c) is
--
--    alpha * a * b + beta * c,
--
matrixMultiply :: Double                          -- alpha
               -> MU.MMatrix VU.MVector s Double  -- matrix a
               -> MU.MMatrix VU.MVector s Double  -- matrix b
               -> Double                          -- beta
               -> MU.MMatrix VU.MVector s Double  -- matrix c, converted to
               -> ST s ()                         -- alpha * (a * b) + beta * c
matrixMultiply alpha a b beta c = do
    let
        (ra, ca) = MU.dim a
        (rb, cb) = MU.dim b
        (rc, cc) = MU.dim c

    -- check compatibility of dimensions
    if ca == rb && rc >= ra && cc >= cb
        then if alpha == 0
                then if beta == 0
                     then numLoop 0 (ra - 1) $ \i ->
                          numLoop 0 (cb - 1) $ \j -> do
                              MU.unsafeWrite c (i, j) 0
                     else numLoop 0 (ra - 1) $ \i ->
                          numLoop 0 (cb - 1) $ \j -> do
                              cij <- MU.unsafeRead c (i, j)
                              MU.unsafeWrite c (i, j) (beta * cij)
                else if beta == 0
                     then numLoop 0 (ra - 1) $ \i ->
                          numLoop 0 (cb - 1) $ \j ->
                          numLoop 0 (ca - 1) $ \k -> do
                              aik <- MU.unsafeRead a (i, k)
                              bkj <- MU.unsafeRead b (k, j)
                              MU.unsafeWrite c (i, j) (alpha * aik * bkj)
                     else numLoop 0 (ra - 1) $ \i ->
                          numLoop 0 (cb - 1) $ \j -> do
                              cij <- MU.unsafeRead c (i, j)
                              numLoop 0 (ca - 1) $ \k -> do
                                  aik <- MU.unsafeRead a (i, k)
                                  bkj <- MU.unsafeRead b (k, j)
                                  MU.unsafeWrite c (i, j) (alpha * aik * bkj + beta * cij)
        else error "incompatible dimensions"


-- | O(1) Extract sub matrix
--
subMatrix :: (Int, Int)  -- ^ upper left corner of the submatrix
          -> (Int, Int)  -- ^ bottom right corner of the submatrix
          -> MU.MMatrix v s a -> MU.MMatrix v s a
{-# INLINE subMatrix #-}
subMatrix (i,j) (i',j') (MU.MMatrix _ n tda offset vec)
    | m' <= 0 || n' <= 0 = error "incorrect dimensions in subMatrix"
    | otherwise = MU.MMatrix m' n' tda offset' vec
  where
    m' = i' - i + 1
    n' = j' - j + 1
    offset' = offset + i * n + j


dotProduct :: (Num a, VU.Unbox a) => VU.MVector s a -> VU.MVector s a -> ST s a
{-# INLINE dotProduct #-}
dotProduct v1 v2 = numLoopState 0 (VU.length v1 - 1) 0 $ \acc i -> do
    v1' <- VU.unsafeRead v1 i
    v2' <- VU.unsafeRead v2 i
    return $ v1' * v2' + acc


testMul :: M.Matrix V.Vector Double
        -> M.Matrix V.Vector Double
        -> M.Matrix V.Vector Double
testMul a b = runST $ do
    a' <- M.thaw a
    b' <- M.thaw b
    c' <- MU.new (M.rows a, M.cols b)
    matrixMultiply 2.0 a' b' 0.0 c'
    M.unsafeFreeze c'


testSub :: M.Matrix V.Vector Double -> M.Matrix V.Vector Double
testSub m = runST $ do
     m' <- M.thaw m
     let ms = subMatrix (1, 1) (1, 1) m'
     m'' <- M.unsafeFreeze ms
     return m''


testSwap :: M.Matrix V.Vector Double -> V.Vector Int -> M.Matrix V.Vector Double
testSwap m p = runST $ do
     m' <- M.thaw m
     p' <- V.thaw p
     let ms = subMatrix (1, 1) (2, 2) m'
     rowSwap ms 0 1 p' Up
     m'' <- M.unsafeFreeze ms
     return m''



data Direction = Up | Down deriving (Eq, Read, Show)

-- rowSwap swaps two rows.  Note that the pivot vector is not
-- arranged as a permutation vector (i.e., the entry at index
-- i corresponding to the row swapped with i), but in NAG pivot
-- format, in which the i-th entry gives the row number that
-- was swapped with i when row i was processed).  An easy way
-- to distinuguish the formats is that in a permutation vector,
-- no entry can be repeated, which in NAG pivot format entries
-- may repeat.
--
-- Note that in either format, an entry which is the same as its
-- index indicates a row that was not swapped.
--
rowSwap :: MU.MMatrix VU.MVector s Double
        -> Int
        -> Int
        -> VU.MVector s Int
        -> Direction
        -> ST s ()
rowSwap a firstPivot lastPivot pivots _ = do
    let
        (nr, nc) = MU.dim a

    if nr /= VU.length pivots
       then error "length of pivot vector must equal number of rows"
       else numLoop firstPivot lastPivot $ \i -> do
            ip <- VU.unsafeRead pivots i
            if ip /= i
               then numLoop 0 (nc - 1) $ \k -> do
                  temp <- MU.unsafeRead a (i,  k)
                  aipk <- MU.unsafeRead a (ip, k)
                  MU.unsafeWrite a (i,  k) aipk
                  MU.unsafeWrite a (ip, k) temp
               else return ()


triangularSolve :: MU.MMatrix VU.MVector s Double
                -> VU.MVector s Double
                -> VU.MVector s Double
                -> ST s ()
triangularSolve a b x = undefined


-- Solve the system of equations Ax = b, given A as a packed LU decomposition
-- and a row permutation vector.
--
-- The arguments are structure so one can solve the linear system ax = b
-- using
--
--    x = luSolve (luFactor a) b
--
luSolve :: (M.Matrix V.Vector Double,    -- matrix A, as a packed LU decomposition
            V.Vector Int)                -- row permutation vector
        -> V.Vector Double               -- vector b
        -> V.Vector Double               -- vector x
luSolve (a, perm) b = runST $ do undefined



-- multStd__ :: Num a => Matrix a -> Matrix a -> Matrix a
-- {-# INLINE multStd__ #-}
-- multStd__ a b = matrix r c $ \(i,j) -> dotProduct (V.unsafeIndex avs $ i - 1) (V.unsafeIndex bvs $ j - 1)
--   where
--     r = nrows a
--     avs = V.generate r $ \i -> getRow (i+1) a
--     c = ncols b
--     bvs = V.generate c $ \i -> getCol (i+1) b

-- dotProduct :: Num a => V.Vector a -> V.Vector a -> a
-- {-# INLINE dotProduct #-}
-- dotProduct v1 v2 = numLoopFold 0 (V.length v1 - 1) 0 $
--   \r i -> V.unsafeIndex v1 i * V.unsafeIndex v2 i + r
