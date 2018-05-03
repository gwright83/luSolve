{-# LANGUAGE AllowAmbiguousTypes #-}
--
-- LUSolve.hs
--
-- Solve an LU factored system of equations.
--

module Numeric.LinearAlgebra.LUSolve where

import           Control.Loop                      (numLoop, numLoopState)
import           Control.Monad.ST                  (ST, runST)
import Control.Monad.ST.Unsafe (unsafeIOToST)
import qualified Data.Matrix.Dense.Generic         as M
import qualified Data.Matrix.Dense.Generic.Mutable as MU
import qualified Data.Vector.Unboxed               as V
import qualified Data.Vector.Unboxed.Mutable       as VU
import Debug.Trace


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
luFactor :: M.Matrix V.Vector Double -> (M.Matrix V.Vector Double, V.Vector Int)
luFactor aOrig = runST $ do
    let
        (m, n) = M.dim aOrig
        mnMin  = min m n

    a      <- M.thaw aOrig
    pivots <- VU.new n

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

    mPeek "a: " a
    a'      <- M.freeze a
    pivots' <- V.freeze pivots

    return (a', pivots')


luFactor_ :: MU.MMatrix VU.MVector s Double
          -> VU.MVector s Int
          -> ST s ()
luFactor_ a pivots = do
    let
        (m, n) = MU.dim a
        n'     = n `div` 2

        mm1  = m  - 1
        nm1  = n  - 1
        npm1 = n' - 1

    mPeek "a: " a
    --mvPeek a
    vPeek pivots

    if n == 1
       then do
            pivotAndScale a pivots
            return ()
       else do
            let

               aLeft  = subMatrix (0, 0)  (mm1, npm1) a
               aRight = subMatrix (0, n') (mm1, nm1)  a

               aTopLeft     = subMatrix (0,  0)  (npm1, npm1) a
               aTopRight    = subMatrix (0,  n') (npm1, nm1)  a
               aBottomLeft  = subMatrix (n', 0)  (mm1,  npm1) a
               aBottomRight = subMatrix (n', n') (mm1,  nm1)  a

               pivotsTop    = VU.slice 0       n'  pivots
               pivotsBottom = VU.slice n' (n - n') pivots

            luFactor_ aLeft  pivotsTop -- was pivotsTop

        --mPeek "aRight: " aRight
        --vPeek pivotsTop

            rowSwap   aRight pivotsTop -- was pivotsTop

            mPeek "aRight': " aRight
            vPeek pivotsTop

            mPeek "aTopRight in: " aTopRight
            triangularSolve aTopLeft  aTopRight
            mPeek "aTopRight out: " aTopRight

            mPeek "A_BR: " aBottomRight
            mPeek "A_BL: " aBottomLeft
            mPeek "A_TR: " aTopRight
            matrixMultiply (-1.0) aBottomLeft aTopRight 1.0 aBottomRight
            mPeek "A_BR': " aBottomRight

            luFactor_ aBottomRight pivotsBottom

        --mPeek "aBottomLeft: " aBottomLeft
        --vPeek pivotsBottom

            rowSwap   aBottomLeft  pivotsBottom

        --mPeek "aBottomLeft': " aBottomLeft
        --vPeek pivotsBottom

            adjustPivots pivotsBottom n'


mPeek :: String -> MU.MMatrix VU.MVector s Double -> ST s ()
mPeek str a = do
    let
        (m, n) = MU.dim a
    numLoop 0 (m - 1) $ \i -> do
      unsafeIOToST (putStr "\n")
      numLoop 0 (n - 1) $ \j -> do
        aij <- MU.read a (i, j)
        unsafeIOToST (putStr (str ++ (show aij) ++ "  "))
    unsafeIOToST (putStr "\n")


mvPeek :: MU.MMatrix VU.MVector s Double -> ST s ()
mvPeek (MU.MMatrix _ _ _ _ v)  = do
    let
        n = VU.length v
    unsafeIOToST (putStr "\n")
    numLoop 0 (n - 1) $ \i -> do
        vi <- VU.read v i
        unsafeIOToST (putStrLn (show vi))


vPeek :: VU.MVector s Int -> ST s ()
vPeek v = do
    let
        n = VU.length v
    unsafeIOToST (putStr "\n")
    numLoop 0 (n - 1) $ \i -> do
        vi <- VU.read v i
        unsafeIOToST (putStrLn (show vi))


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
    v1' <- VU.read v1 i
    v2' <- VU.read v2 i
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
     rowSwap ms p'
     m'' <- M.unsafeFreeze ms
     return m''


testPivotAndScale :: M.Matrix V.Vector Double
                  -> V.Vector Int
                  -> (M.Matrix V.Vector Double, V.Vector Int)
testPivotAndScale m p = runST $ do
     m' <- M.thaw m
     p' <- V.thaw p
     pivotAndScale m' p'
     m'' <- M.unsafeFreeze m'
     p'' <- V.unsafeFreeze p'
     return (m'', p'')


testTriangularSolve :: M.Matrix V.Vector Double
                    -> M.Matrix V.Vector Double
                    -> M.Matrix V.Vector Double
testTriangularSolve a b = runST $ do
     a' <- M.thaw a
     b' <- M.thaw b
     triangularSolve a' b'
     b'' <- M.unsafeFreeze b'
     return b''


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
        -> VU.MVector s Int
        -> ST s ()
rowSwap a pivots = do
    let
        (_, nc) = MU.dim a
        nPivots = VU.length pivots

    numLoop 0 (nPivots - 1) $ \i -> do
        ip <- VU.read pivots i
        if ip /= i
           then numLoop 0 (nc - 1) $ \k -> do
              let
                  i'  = i
                  ip' = ip
              temp <- MU.read a (i',  k)
              aipk <- MU.read a (ip', k)
              MU.write a (i',  k) aipk
              MU.write a (ip', k) temp
           else return ()


pivotAndScale :: MU.MMatrix VU.MVector s Double
              -> VU.MVector s Int
              -> ST s ()
pivotAndScale a pivots = do
    ip   <- findPivot a
    temp <- MU.read a (0,  0)
    aip  <- MU.read a (ip, 0)
    MU.write a (0,  0) aip
    MU.write a (ip, 0) temp
    VU.write pivots 0 ip

    -- Scale the elememts below the first.
    let
        (nr, _) = MU.dim a

    numLoop 1 (nr - 1) $ \k -> do
      ak <- MU.read a (k, 0)
      MU.write a (k, 0) (ak / aip)


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
        ip <- VU.read pivots i
        VU.write pivots i (ip + offset)


-- TriangularSolve solves the linear system AX = B when A is upper
-- triangular.  The matrix B is overwritten, column by column, by
-- the solution matrix X.
--
triangularSolve :: MU.MMatrix VU.MVector s Double
                -> MU.MMatrix VU.MVector s Double
                -> ST s ()
triangularSolve a b = do
    let
        (_, m) = MU.dim a
        (_, n) = MU.dim b

    numLoop 0 (n - 1) $ \j ->
      numLoop 0 (m - 1) $ \k -> do
        bkj <- MU.read b (k, j)
        if bkj == 0
           then return ()
           else numLoop (k + 1) (m - 1) $ \i -> do
            bij <- MU.read b (i, j)
            aik <- MU.read a (i, k)
            MU.write b (i, j) (bij - bkj * aik)


-- Should be correct, but still needs to be tested.
--
findPivot :: MU.MMatrix VU.MVector s Double -> ST s Int
findPivot m = do
    let
        (nr, _) = MU.dim m

        go aval k idx = do
        if k >= nr
           then return idx
           else do
            v <- MU.read m (k, 0)
            if aval < abs v
               then go (abs v) (k + 1) k
               else go  aval   (k + 1) idx

    val <- MU.read m (0, 0)
    piv <- go (abs val) 1 0
    return piv


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


testMat = M.fromLists [[0.772386, 0.499327, 0.189312],
                       [0.759731, 0.799350, 0.682719],
                       [0.456574, 0.636521, 0.003035],
                       [0.014020, 0.636044, 0.990054]]

test :: (M.Matrix V.Vector Double, V.Vector Int)
test = luFactor testMat
