--
-- Spec file for luSolve.hs
--
-- Tests the exposed LU decomposition and linear system solver.
--
--
module Numeric.LinearAlgebra.LUSolveSpec where

import           Test.Hspec

import           Data.Matrix.Generic           as M hiding (zipWith)
import qualified Data.Vector.Unboxed           as V
import           Numeric.LinearAlgebra.LUSolve
import           System.Random


-- A useful, but strangely non-standard function: given an integer n
-- and a list, break the input list up into lists of length n and
-- return the list of those lists.
--
-- The name reflects that it "bundles" a list into lists of length n.
--
-- If the length of the input list is not a multiple of n, the final
-- sublist will be shorter than n elements.
--
-- The definition here is works for infinite lists, which is perhaps
-- the most common use-case for this function.
--
bundle :: Int -> [ a ] -> [[ a ]]
bundle _ [] = []
bundle n xs = take n xs : bundle n (drop n xs)


mVals :: [ Double ]
mVals = randoms (mkStdGen 1)


-- Make an infinite list of random n * n matrices.
--
randomSquareMatrices :: Int -> [ M.Matrix V.Vector Double ]
randomSquareMatrices n = Prelude.map (\vs -> M.fromLists (bundle n vs)) (bundle (n * n) mVals)

vVals :: [ Double ]
vVals = randoms (mkStdGen 2)

-- Make an infinite list of random n-vectors.
--
randomColumnVectors :: Int -> [ M.Matrix V.Vector Double ]
randomColumnVectors n = Prelude.map (\vs -> M.fromLists (bundle 1 vs )) (bundle n vVals)


-- Solve the linear system Ax = b for x, then compute the maximum
-- absolute value of the vector of differences Ax - b.  The later
-- is the value of the function.
--
checkLinearSystemSolution :: M.Matrix V.Vector Double  -- matrix A
                          -> M.Matrix V.Vector Double  -- right hand vector b
                          -> Double                    -- max (abs (Ax - b))
checkLinearSystemSolution a b = let
    x  = luSolve (luFactor a) b
    b' = matrixMultiply a x
    in
      maxAbsDiff b b'


-- Multiply two matrices
--
matrixMultiply :: M.Matrix V.Vector Double
               -> M.Matrix V.Vector Double
               -> M.Matrix V.Vector Double
matrixMultiply a b = let
    (ra, ca) = dim a
    (rb, cb) = dim b
    in
      if ca /= rb
      then error "incompatible dimensions in matrixMultiply"
      else fromLists [[ sum [ a ! (i, k) * b ! (k, j) | k <- [0 .. (ca - 1)]]
                                                      | j <- [0 .. (cb - 1)]]
                                                      | i <- [0 .. (ra - 1)]]


-- Return the largest difference (by absolute value) between
-- two matrices.
--
maxAbsDiff :: M.Matrix V.Vector Double
           -> M.Matrix V.Vector Double
           -> Double
maxAbsDiff (M.Matrix m n _ _ v) (M.Matrix m' n' _ _ v') =
    if m == m' && n == n'
    then V.maximum $ V.zipWith (\y y' -> abs (y - y')) v v'
    else error "unequal matrix dimensions in maxAbsDiff"


-- Return the largest difference (by relative value) between
-- two matrices.
--
maxRelDiff :: M.Matrix V.Vector Double
           -> M.Matrix V.Vector Double
           -> Double
maxRelDiff (Matrix m n _ _ v) (Matrix m' n' _ _ v') =
    if m == m' && n == n'
    then V.maximum $ V.zipWith (\y y' -> abs (2 * (y - y') / (abs y + abs y'))) v v'
    else error "unequal matrix dimensions in maxRelDiff"


-- Threshold for accepting max (abs (Ax - b)) as correct.
--
matrixSpecAbsDiff :: Double
matrixSpecAbsDiff = 1.0e-9


-- Generate n linear systems of size 'size' and solve them by luSolve.
-- Return a list of the maximum absolute value of the error (max (abs (Ax - b)))
-- for each case.
--
runCheck :: Int         -- number of cases to run
         -> Int         -- number of equations in each linear system
         -> [ Double ]  -- list of maximum absolute values of the errors
runCheck n size = take n $ zipWith checkLinearSystemSolution (randomSquareMatrices size)
                                                             (randomColumnVectors  size)


-- The tests themselves
--
spec :: Spec
spec = do
    -- describe "Numeric.LinearAlgebra.LUSolve.luSolve" $ do
    --     it "Solve 1000 random 10 * 10 systems" $ do
    --         (maximum $ runCheck 1000 10)
    --             `shouldSatisfy` (< matrixSpecAbsDiff)

    -- describe "Numeric.LinearAlgebra.LUSolve.luSolve" $ do
    --     it "Solve 100 random 50 * 50 systems" $ do
    --         (maximum $ runCheck 100 50)
    --             `shouldSatisfy` (< matrixSpecAbsDiff)

    -- describe "Numeric.LinearAlgebra.LUSolve.luSolve" $ do
    --     it "Solve 100 random 100 * 100 systems" $ do
    --         (maximum $ runCheck 100 100)
    --             `shouldSatisfy` (< matrixSpecAbsDiff)

    -- describe "Numeric.LinearAlgebra.LUSolve.luSolve" $ do
    --     it "Solve 100 random 500 * 500 systems" $ do
    --         (maximum $ runCheck 100 500)
    --             `shouldSatisfy` (< matrixSpecAbsDiff)

    describe "Numeric.LinearAlgebra.LUSolve.luSolve" $ do
        it "Solve 100 random 1000 * 1000 systems" $ do
            (maximum $ runCheck 100 1000)
                `shouldSatisfy` (< matrixSpecAbsDiff)


main :: IO ()
main = hspec spec
