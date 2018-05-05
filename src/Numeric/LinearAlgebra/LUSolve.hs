-- |
-- Module      : Numeric.LinearAlgebra.LUSolve
-- Copyright   : (c) Gregory Wright 2018
-- License     : BSD-style
--
-- Maintainer  : Gregory Wright <gwright@antiope.com.
-- Stability   : experimental
-- Portability : non-portable
--
-- The LUSolve library implements Crout's algorithm for in-place
-- LU decomposition in pure Haskell.  A function using the
-- LU-factored matrix to solve systems of linear equations is
-- also provided.
--
-- The version of Crout's algorithm is that given by F.G. Gustavson,
-- /IBM Journal of Research and Development/, Vol. 41, No. 6, 1997,
-- pp. 737-755.  This is a recursive, in-place, procedure.  Here
-- it is implemented using a mutable matrix in the ST monad.
-- Partial (row) pivoting is used for numerical stability.
--

module Numeric.LinearAlgebra.LUSolve (
    -- * LU Decomposition
    luFactor,
    luFactor_,
    -- * Solving Linear Systems
    luSolve
    ) where

import           Control.Loop                      (numLoop)
import           Control.Monad.ST                  (ST, runST)
import qualified Data.Matrix.Dense.Generic         as M
import qualified Data.Matrix.Dense.Generic.Mutable as MU
import qualified Data.Vector.Unboxed               as V
import qualified Data.Vector.Unboxed.Mutable       as VU


-- | LU Decomposition
--
-- Factor a rectangular matrix A into PA = LU, where L is lower triangular,
-- U is upper triangular and P is the permuation matrix (represented as
-- a vector containing the location of the nonzero column) that describes
-- the row interchanges (partial pivoting).
--
-- To maintain the expected Haskell API, i.e., an immutable input matrix,
-- the first thing to do is to copy the original immutable
-- matrix into a mutable one.  The LU factorization runs efficiently in place
-- using a mutable matrix. At the end of the calculation, the
-- matrix will be frozen (i.e., made immutable again).
--
-- The factorization takes place in eiher two steps or a single step.
-- If the number of rows of aOrig (m) is less than the number of columns (n),
-- the m * m square matrix is factored first, then the remaining m * (n - m) piece.
-- The the number of rows is greater than or equal to the number of columns,
-- one invocation of luFactor_ is all that is needed.
--
-- The LU factored matrix is returned in packed format. The upper triangular
-- part is the U matrix.  The lower triangular part beneath the main diagonal
-- is the L matrix without its diagonal entries, which are omitted since they
-- are known to be 1. This is the traditional way of storing the LU decomposition
-- and luSolve understand this format.
--
-- The returned pivot vector is not a permutaion vector, but instead is
-- in "NAG pivot format".  Reading the vector from top to bottom
-- (equivalently, frmo left to right), the current entry specifies which
-- row to swap with the current row.  Note that unlike in a permutation vector,
-- where each entry is the index of the nonzero entry of the permutation
-- matrix, a NAG format pivot vector can have repeated entries.
--
luFactor ::  M.Matrix V.Vector Double   -- ^ matrix A
         -> (M.Matrix V.Vector Double,  -- ^ LU decomposition of A
             V.Vector Int)              -- ^ row pivots
luFactor aOrig = runST $ do
    let
        (m, n) = M.dim aOrig
        mnMin  = min m n

    a      <- M.thaw aOrig
    pivots <- VU.unsafeNew mnMin

    let
        a' = subMatrix (0, 0) (m - 1, mnMin - 1) a

    luFactor_ a' pivots

    if m >= n
       then do return ()
       else do
        let
            aLeft  = subMatrix (0, 0) (m - 1, m - 1) a
            aRight = subMatrix (0, m) (m - 1, n - 1) a

        rowSwap aRight pivots
        triangularSolve aLeft aRight

    aFactored <- M.freeze a
    pivots'   <- V.freeze pivots

    return (aFactored, pivots')


-- | The luFactor function takes a mutable matrix and replaces
-- it with its LU decomposition.  An unitialized pivot vector
-- is replaced with the row pivots, in NAG pivot format.
--
luFactor_ :: MU.MMatrix VU.MVector s Double  -- ^ matrix A, overwritten by LU
          -> VU.MVector s Int                -- ^ row pivots
          -> ST s ()
luFactor_ a pivots = do
    let
        (m, n) = MU.dim a
        n'     = n `div` 2

    if n == 1
       then pivotAndScale a pivots
       else do
            let

               aLeft  = subMatrix (0, 0)  (m - 1, n' - 1) a
               aRight = subMatrix (0, n') (m - 1, n  - 1) a

               aTopLeft     = subMatrix (0,  0)  (n' - 1, n' - 1) a
               aTopRight    = subMatrix (0,  n') (n' - 1, n  - 1) a
               aBottomLeft  = subMatrix (n', 0)  (m  - 1, n' - 1) a
               aBottomRight = subMatrix (n', n') (m  - 1, n  - 1) a

               pivotsTop    = VU.slice 0       n'  pivots
               pivotsBottom = VU.slice n' (n - n') pivots

            luFactor_ aLeft  pivotsTop
            rowSwap   aRight pivotsTop
            triangularSolve aTopLeft  aTopRight
            matrixMultiply (-1.0) aBottomLeft aTopRight 1.0 aBottomRight
            luFactor_ aBottomRight pivotsBottom
            rowSwap   aBottomLeft  pivotsBottom

            -- Add an offset to pivotsBottom it entries refer to the
            -- row number of the original matrix, rather than the
            -- submatrix.
            adjustPivots pivotsBottom n'


-- |  Solve the system of equations Ax = b, given A as a packed LU decomposition
-- and a row permutation vector in NAG pivot format.
--
-- The arguments are structured so one can solve the linear system Ax = b
-- using
--
--    x = luSolve (luFactor a) b
--
luSolve :: (M.Matrix V.Vector Double,    -- ^ matrix A, as a packed LU decomposition
            V.Vector Int)                -- ^ row pivots
        -> V.Vector Double               -- ^ vector b
        -> V.Vector Double               -- ^ vector x
luSolve (_, _) _ = runST $ do undefined



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
                              MU.write c (i, j) 0
                     else numLoop 0 (ra - 1) $ \i ->
                          numLoop 0 (cb - 1) $ \j -> do
                              cij <- MU.read c (i, j)
                              MU.write c (i, j) (beta * cij)
                else if beta == 0
                     then numLoop 0 (ra - 1) $ \i ->
                          numLoop 0 (cb - 1) $ \j ->
                          numLoop 0 (ca - 1) $ \k -> do
                              aik <- MU.read a (i, k)
                              bkj <- MU.read b (k, j)
                              MU.write c (i, j) (alpha * aik * bkj)
                     else numLoop 0 (ra - 1) $ \i ->
                          numLoop 0 (cb - 1) $ \j -> do
                              cij <- MU.read c (i, j)
                              numLoop 0 (ca - 1) $ \k -> do
                                  aik <- MU.read a (i, k)
                                  bkj <- MU.read b (k, j)
                                  MU.write c (i, j) (alpha * aik * bkj + beta * cij)
        else error "incompatible dimensions"


-- Extract sub matrix
--
-- This is the same as the subMatrix function exported by the
-- matrices library, but for a mutable matrix.
--
subMatrix :: (Int, Int)  -- ^ upper left corner of the submatrix
          -> (Int, Int)  -- ^ bottom right corner of the submatrix
          -> MU.MMatrix v s a -> MU.MMatrix v s a
{-# INLINE subMatrix #-}
subMatrix (i,j) (i',j') (MU.MMatrix _ _ tda offset vec)
    | m' <= 0 || n' <= 0 = error "incorrect dimensions in subMatrix"
    | otherwise = MU.MMatrix m' n' tda offset' vec
  where
    m' = i' - i + 1
    n' = j' - j + 1
    offset' = offset + i * tda + j


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
-- index indicates a row that is not swapped.
--
rowSwap :: MU.MMatrix VU.MVector s Double
        -> VU.MVector s Int
        -> ST s ()
rowSwap a pivots = do
    let
        (_, nc) = MU.dim a
        nPivots = VU.length pivots

    numLoop 0 (nPivots - 1) $ \i -> do
        ip <- VU.unsafeRead pivots i
        if ip /= i
           then numLoop 0 (nc - 1) $ \k -> do
              let
                  i'  = i
                  ip' = ip
              temp <- MU.unsafeRead a (i',  k)
              aipk <- MU.unsafeRead a (ip', k)
              MU.unsafeWrite a (i',  k) aipk
              MU.unsafeWrite a (ip', k) temp
           else return ()


-- pivotAndScale computes the LU decompostion of a matrix
-- with a single column.
--
pivotAndScale :: MU.MMatrix VU.MVector s Double
              -> VU.MVector s Int
              -> ST s ()
pivotAndScale a pivots = do
    ip   <- findPivot a
    temp <- MU.unsafeRead a (0,  0)
    aip  <- MU.unsafeRead a (ip, 0)
    MU.unsafeWrite a (0,  0) aip
    MU.unsafeWrite a (ip, 0) temp
    VU.unsafeWrite pivots 0 ip

    -- Scale the elememts below the first.
    let
        (nr, _) = MU.dim a

    numLoop 1 (nr - 1) $ \k -> do
      ak <- MU.unsafeRead a (k, 0)
      MU.unsafeWrite a (k, 0) (ak / aip)


-- Given a pivot vector, add a constant to each element.
-- This is used to shift the pivot vector entries from referring
-- to the local submatrix to the global matrix.
--
adjustPivots :: VU.MVector s Int
             -> Int
             -> ST s ()
adjustPivots pivots offset = do
    let
        nPivots = VU.length pivots

    numLoop 0 (nPivots - 1) $ \i -> do
        ip <- VU.unsafeRead pivots i
        VU.unsafeWrite pivots i (ip + offset)


-- TriangularSolve solves the linear system AX = B where A is upper
-- triangular.  The matrix B is overwritten, column by column, by
-- the solution matrix X.
--
triangularSolve :: MU.MMatrix VU.MVector s Double   -- matrix A
                -> MU.MMatrix VU.MVector s Double   -- matrix B
                -> ST s ()
triangularSolve a b = do
    let
        (_, m) = MU.dim a
        (_, n) = MU.dim b

    numLoop 0 (n - 1) $ \j ->
      numLoop 0 (m - 1) $ \k -> do
        bkj <- MU.unsafeRead b (k, j)
        if bkj == 0
           then return ()
           else numLoop (k + 1) (m - 1) $ \i -> do
            bij <- MU.unsafeRead b (i, j)
            aik <- MU.unsafeRead a (i, k)
            MU.unsafeWrite b (i, j) (bij - bkj * aik)


-- Return the index of the matrix element with the largest absolute
-- value in the first column.
--
findPivot :: MU.MMatrix VU.MVector s Double -> ST s Int
findPivot m = do
    let
        (nr, _) = MU.dim m

        go aval k idx = do
            if k >= nr
                then return idx
                else do
                    v <- MU.unsafeRead m (k, 0)
                    if aval < abs v
                        then go (abs v) (k + 1) k
                        else go  aval   (k + 1) idx

    val <- MU.unsafeRead m (0, 0)
    piv <- go (abs val) 1 0
    return piv


testMat :: M.Matrix V.Vector Double
testMat = M.fromLists [[0.772386, 0.499327, 0.189312],
                       [0.759731, 0.799350, 0.682719],
                       [0.456574, 0.636521, 0.003035],
                       [0.014020, 0.636044, 0.990054]]

testMat' :: M.Matrix V.Vector Double
testMat' = M.fromLists [[0.772386, 0.499327, 0.189312, 0.014020],
                        [0.759731, 0.799350, 0.682719, 0.636044],
                        [0.456574, 0.636521, 0.003035, 0.990054]]


_test :: (M.Matrix V.Vector Double, V.Vector Int)
_test = luFactor testMat

_test' :: (M.Matrix V.Vector Double, V.Vector Int)
_test' = luFactor testMat'
