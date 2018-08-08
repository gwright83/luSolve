-- |
-- Module      : Numeric.LinearAlgebra.LUSolve
-- Copyright   : (c) Gregory Wright 2018
-- License     : BSD-style
--
-- Maintainer  : Gregory Wright <gwright@antiope.com>
-- Stability   : experimental
-- Portability : non-portable
--
-- The LUSolve library implements Crout's algorithm for in-place
-- LU decomposition in pure Haskell.  A function using the
-- LU-factored matrix to solve a system of linear equations is
-- also provided.
--
-- The version of Crout's algorithm is that given by F.G. Gustavson,
-- /IBM Journal of Research and Development/, Vol. 41, No. 6, 1997,
-- pp. 737-755.  It is a recursive, in-place, procedure.  This
-- Haskell implementation uses a mutable matrix in the ST monad.
-- Partial (row) pivoting provides numerical stability.
--

module Numeric.LinearAlgebra.LUSolve (
    -- * LU Decomposition
    luFactor,
    luFactor_,
    -- * Solving Linear Systems
    luSolve
    ) where

import           Control.Loop                      (forLoop, numLoop,
                                                    numLoopState)
import           Control.Monad                     (when)
import           Control.Monad.ST                  (ST, runST)
import qualified Data.Matrix.Dense.Generic         as M
import qualified Data.Matrix.Dense.Generic.Mutable as MU
import           Data.STRef.Strict
import qualified Data.Vector.Unboxed               as V
import qualified Data.Vector.Unboxed.Mutable       as VU


-- | LU Decomposition
--
-- Factor a rectangular matrix A into PA = LU, where L is lower triangular,
-- U is upper triangular and P is the permuation matrix that gives
-- the row interchanges (partial pivoting).  Because the permutation matrix
-- is sparse, it is stored in a special format, described below.
--
-- To maintain the expected Haskell API, i.e., an immutable input matrix,
-- the original matrix is copied to a mutable one. LU factorization runs
-- efficiently in place using a mutable matrix. At the end, the mutable
-- matrix is frozen (made immutable).
--
-- The factorization takes place in either two steps or one. If the number
-- of rows of A, m, is less than the number of columns, n, the m * m square
-- matrix is factored first, followed by the remaining m * (n - m) piece.
-- If the the number of rows is greater than or equal to the number of columns,
-- only a single invocation of luFactor_ is required.
--
-- The LU factored matrix is returned in a packed format. The upper triangular
-- part is the U matrix.  The lower triangular part beneath the main diagonal
-- is the L matrix without its diagonal entries, which are omitted since they
-- are known to be 1. This is the traditional way of storing the LU decomposition
-- and linear system solver luSolve takes this format as input.
--
-- The returned pivot vector is not a permutation vector, but instead is
-- in "NAG pivot format".  Considering the pivot vector as a column and reading
-- sequentially from top to bottom, the current entry specifies which
-- row to swap with the current row.  Note that unlike a permutation vector,
-- in which each element is the (unique) index of a nonzero entry in the
-- permutation matrix, a NAG pivot vector can have repeated entries.
--
-- The parity, equal to (-1)^(number of row interchanges), is also returned.
-- The determinant of a square input matrix is the product of the diagonal
-- entries of its packed LU decomposition and the parity.
--
luFactor ::  M.Matrix V.Vector Double   -- ^ Matrix A
         -> (M.Matrix V.Vector Double,
             V.Vector Int,
             Int)                       -- ^ (LU Decomposition of A,
                                        --    row pivots,
                                        --    parity = (-1)^(number of row interchanges))
luFactor aOrig = runST $ do
    let
        (m, n) = M.dim aOrig
        mnMin  = min m n

    a      <- M.thaw aOrig         -- thaw forces a copy, since we might use aOrig again.
    pivots <- VU.unsafeNew mnMin   -- unsafe, since it not used outside this function.
    parity <- newSTRef 1

    let
        a' = subMatrix (0, 0) (m - 1, mnMin - 1) a

    luFactor_ a' pivots parity

    when (m < n) $ do
        let
            aLeft  = subMatrix (0, 0) (m - 1, m - 1) a
            aRight = subMatrix (0, m) (m - 1, n - 1) a

        rowSwap aRight pivots
        triangularSolve Lower Unit aLeft aRight

    aFactored <- M.unsafeFreeze a
    pivots'   <- V.unsafeFreeze pivots
    parity'   <- readSTRef parity

    return (aFactored, pivots', parity')


-- | The luFactor_ function takes a mutable matrix and replaces
-- it with its LU decomposition.  An unitialized pivot vector
-- is replaced with the row pivots, in NAG pivot format.
--
luFactor_ :: MU.MMatrix VU.MVector s Double  -- ^ matrix A, overwritten by LU
          -> VU.MVector s Int                -- ^ row pivots
          -> STRef s Int                     -- ^ parity
          -> ST s ()
luFactor_ a pivots parity = do
    let
        (m, n) = MU.dim a
        n'     = n `div` 2

    if n == 1
       then pivotAndScale a pivots parity
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

            luFactor_ aLeft  pivotsTop parity
            rowSwap   aRight pivotsTop
            triangularSolve Lower Unit aTopLeft  aTopRight
            matrixMultiply (-1.0) aBottomLeft aTopRight 1.0 aBottomRight
            luFactor_ aBottomRight pivotsBottom parity
            rowSwap   aBottomLeft  pivotsBottom

            -- Add an offset to pivotsBottom it entries refer to the
            -- row number of the original matrix, rather than the
            -- submatrix.
            adjustPivots pivotsBottom n'


-- |  Solve the system of equations AX = B, given A as a packed LU decomposition
-- of a  square matrix and a row permutation vector in NAG pivot format.
-- Note that X and B are matrices, not vectors.  This allows solving for multiple
-- right hand sides simultanously.
--
-- The arguments are structured so one can solve the linear system AX = B
-- using
--
-- @
--   x = luSolve (luFactor a) b
-- @
luSolve :: (M.Matrix V.Vector Double,
            V.Vector Int,
            Int)                         -- ^ (Matrix A as a packed LU decompostion,
                                         --    row pivots,
                                         --    parity (ignored here))
        -> M.Matrix V.Vector Double      -- ^ matrix B
        -> M.Matrix V.Vector Double      -- ^ matrix X
luSolve (lu, pivots, _) b = runST $ do
    let
        (m, n) = M.dim lu
        (l, _) = M.dim b

    x <- M.thaw b

    if m /= n
       then error "under- (or over-) determined system in luSolve"
       else if n /= l
        then error "incompatible dimensions in luSolve_"
        else luSolve_ lu pivots x

    x' <- M.unsafeFreeze x
    return x'


-- | luSolve_ does the work of solving the system of linear equations Ax = b.
-- The input matrix b is overwritten by the solution x.  Dimensions
-- are not checked; that should be done by the calling function.
--
luSolve_ :: M.Matrix V.Vector Double        -- ^ matrix A, as a packed LU decomposition
         -> V.Vector Int                    -- ^ pivot vector
         -> MU.MMatrix VU.MVector s Double  -- ^ right hand side b, overwritten by x
         -> ST s ()
luSolve_ lu pivots b = do
    lu'     <- M.unsafeThaw lu       -- Can be unsafe, since it is only read.
    pivots' <- V.unsafeThaw pivots   -- Ditto.
    rowSwap b pivots'
    triangularSolve Lower Unit    lu' b
    triangularSolve Upper NonUnit lu' b


-- This is a generic mutable matrix multiply.  The mutable references
-- need not be distinct, allowing it to be used as part of an in-place
-- algorithm like LU factorization.
--
-- What is computed (and stored in the matrix c) is
--
--    alpha * a * b + beta * c
--
-- so the operation is a matrix 'fused multiply add'.
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
    if ca == rb && rc == ra && cc == cb
        then if alpha == 0
                then if beta == 0
                     then numLoop 0 (ra - 1) $ \i ->
                          numLoop 0 (cb - 1) $ \j -> do
                              MU.unsafeWrite c (i, j) 0
                     else numLoop 0 (ra - 1) $ \i ->
                          numLoop 0 (cb - 1) $ \j -> do
                              cij <- MU.unsafeRead c (i, j)
                              MU.unsafeWrite c (i, j) (beta * cij)
                else do
                     if beta == 0
                     then numLoop 0 (ra - 1) $ \i ->
                          numLoop 0 (cb - 1) $ \j -> do
                              s <- numLoopState 0 (ca - 1) 0 $ \s k -> do
                                  aik <- MU.unsafeRead a (i, k)
                                  bkj <- MU.unsafeRead b (k, j)
                                  return $! s + aik * bkj
                              MU.unsafeWrite c (i, j) (alpha * s)
                     else numLoop 0 (ra - 1) $ \i ->
                          numLoop 0 (cb - 1) $ \j -> do
                              s <- numLoopState 0 (ca - 1) 0 $ \s k -> do
                                  aik <- MU.unsafeRead a (i, k)
                                  bkj <- MU.unsafeRead b (k, j)
                                  return $! s + aik * bkj
                              cij <- MU.unsafeRead c (i, j)
                              MU.unsafeWrite c (i, j) (alpha * s + beta * cij)
         else error "incompatible dimensions"

_testMul :: Double
         -> M.Matrix V.Vector Double
         -> M.Matrix V.Vector Double
         -> Double
         -> M.Matrix V.Vector Double
         -> M.Matrix V.Vector Double
_testMul alpha a b beta c = runST $ do
    a' <- M.thaw a
    b' <- M.thaw b
    c' <- M.thaw c
    matrixMultiply alpha a' b' beta c'
    c'' <- M.freeze c'
    return c''

-- Extract a sub matrix
--
-- This is the same as the subMatrix function exported by the
-- matrices library, but for a mutable matrix.
--
subMatrix :: (Int, Int)  -- ^ upper left corner of the submatrix
          -> (Int, Int)  -- ^ bottom right corner of the submatrix
          -> MU.MMatrix v s a
          -> MU.MMatrix v s a
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
-- i corresponding to the row swapped with i), but is in NAG pivot
-- format, in which the i-th entry gives the row number that
-- was swapped with i when row i was processed).  An easy way
-- to distinuguish the formats is that in a permutation vector,
-- every entry must be unique, which in NAG pivot format entries
-- may be duplicated.
--
-- Note that in either format, an entry which is the same as its
-- index indicates a row that is not swapped.
--
rowSwap :: MU.MMatrix VU.MVector s Double
        -> VU.MVector s Int
        -> ST s ()
{-# INLINE rowSwap #-}
rowSwap a pivots = do
    let
        (_, nc) = MU.dim a
        nPivots = VU.length pivots

    numLoop 0 (nPivots - 1) $ \i -> do
        ip <- VU.unsafeRead pivots i
        when (ip /= i) $
           numLoop 0 (nc - 1) $ \k -> do
              let
                  i'  = i
                  ip' = ip
              temp <- MU.unsafeRead a (i',  k)
              aipk <- MU.unsafeRead a (ip', k)
              MU.unsafeWrite a (i',  k) aipk
              MU.unsafeWrite a (ip', k) temp


-- pivotAndScale computes the LU decompostion of a matrix
-- with a single column.
--
pivotAndScale :: MU.MMatrix VU.MVector s Double
              -> VU.MVector s Int
              -> STRef s Int
              -> ST s ()
{-# INLINE pivotAndScale #-}
pivotAndScale a pivots parity = do
    ip   <- findPivot a
    temp <- MU.unsafeRead a (0,  0)
    aip  <- MU.unsafeRead a (ip, 0)

    if aip == 0.0
       then error "zero pivot in pivotAndScale"
       else do
        MU.unsafeWrite a (0,  0) aip
        MU.unsafeWrite a (ip, 0) temp
        VU.unsafeWrite pivots 0  ip

        when (ip /= 0) $ modifySTRef' parity (* (-1))

        -- Scale the elements below the first.
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
{-# INLINE adjustPivots #-}
adjustPivots pivots offset = do
    let
        nPivots = VU.length pivots

    numLoop 0 (nPivots - 1) $ \i -> do
        ip <- VU.unsafeRead pivots i
        VU.unsafeWrite pivots i (ip + offset)


data Triangle = Upper | Lower
    deriving (Eq, Show)

data DiagonalUnit = Unit | NonUnit
    deriving (Eq, Show)


-- TriangularSolve solves the linear system AX = B where A is upper
-- triangular.  The matrix B is overwritten, column by column, by
-- the solution matrix X.
--
triangularSolve :: Triangle                         -- is A upper or lower triangular?
                -> DiagonalUnit                     -- diagonal entries are 1 or not?
                -> MU.MMatrix VU.MVector s Double   -- matrix A
                -> MU.MMatrix VU.MVector s Double   -- matrix B
                -> ST s ()
triangularSolve Lower unit a b = do
    let
        (_, m) = MU.dim a
        (_, n) = MU.dim b

    numLoop 0 (n - 1) $ \j ->
      numLoop 0 (m - 1) $ \k -> do
        bkj <- MU.unsafeRead b (k, j)
        when (bkj /= 0) $ do
            when (unit == NonUnit) $ do
                akk <- MU.unsafeRead a (k, k)
                MU.unsafeWrite b (k, j) (bkj / akk)
            numLoop (k + 1) (m - 1) $ \i -> do
                bij  <- MU.unsafeRead b (i, j)
                aik  <- MU.unsafeRead a (i, k)
                bkj' <- MU.unsafeRead b (k, j)
                MU.unsafeWrite b (i, j) (bij - bkj' * aik)

triangularSolve Upper unit a b = do
    let
        (_, m) = MU.dim a
        (_, n) = MU.dim b

    numLoop 0 (n - 1) $ \j ->
      forLoop (m - 1) (>= 0) (subtract 1) $ \k -> do
        bkj <- MU.unsafeRead b (k, j)
        when (bkj /= 0) $ do
            when (unit == NonUnit) $ do
                akk <- MU.unsafeRead a (k, k)
                MU.unsafeWrite b (k, j)  (bkj / akk)
            numLoop 0 (k - 1) $ \i -> do
                bij  <- MU.unsafeRead b (i, j)
                aik  <- MU.unsafeRead a (i, k)
                bkj' <- MU.unsafeRead b (k, j)
                MU.unsafeWrite b (i, j) (bij - bkj' * aik)


-- Return the index of the matrix element with the largest absolute
-- value in the first column.
--
findPivot :: MU.MMatrix VU.MVector s Double -> ST s Int
{-# INLINE findPivot #-}
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



-- Some test matrices and driver functions for sanity checks
-- during developent.
--
_testMat :: M.Matrix V.Vector Double
_testMat = M.fromLists [[0.772386, 0.499327, 0.189312],
                        [0.759731, 0.799350, 0.682719],
                        [0.456574, 0.636521, 0.003035],
                        [0.014020, 0.636044, 0.990054]]

_testMat' :: M.Matrix V.Vector Double
_testMat' = M.fromLists [[0.772386, 0.499327, 0.189312, 0.014020],
                         [0.759731, 0.799350, 0.682719, 0.636044],
                         [0.456574, 0.636521, 0.003035, 0.990054]]


_testMat'' :: M.Matrix V.Vector Double
_testMat'' = M.fromLists [[4.0, 3.0, 1.0, 3.0],
                          [3.0, 5.0, 6.0, 8.0],
                          [1.0, 6.0, 2.0, 7.0],
                          [3.0, 8.0, 7.0, 9.0]]

_testMat3 :: M.Matrix V.Vector Double
_testMat3 = M.fromLists [[1.0, 2.0], [3.0, 4.0]]

_testMat4 :: M.Matrix V.Vector Double
_testMat4 = M.fromLists [[5.0], [6.0]]

_testMat5 :: M.Matrix V.Vector Double
_testMat5 = M.fromLists [[2.0], [1.0]]

_test :: (M.Matrix V.Vector Double, V.Vector Int, Int)
_test = luFactor _testMat

_test' :: (M.Matrix V.Vector Double, V.Vector Int, Int)
_test' = luFactor _testMat'
