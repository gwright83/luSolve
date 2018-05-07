--
-- Spec file for luSolve.hs
--
-- Tests the exposed LU decomposition and linear system solver.
--
--
module Numeric.LinearAlgebra.LUSolveSpec where

import           Test.Hspec

import           Data.Matrix.Dense.Generic     as M
import qualified Data.Vector.Unboxed           as V
import           Numeric.LinearAlgebra.LUSolve
import           System.Random


bundle :: Int -> [ a ] -> [[ a ]]
bundle _ [] = []
bundle n xs = take n xs : bundle n (drop n xs)

mVals :: [ Double ]
mVals = randoms (mkStdGen 1)


randomSquareMatrices :: Int -> [ M.Matrix V.Vector Double ]
randomSquareMatrices n = Prelude.map (\vs -> M.fromLists (bundle n vs)) (bundle (n * n) mVals)

vVals :: [ Double ]
vVals = randoms (mkStdGen 2)

randomColumnVectors :: Int -> [ M.Matrix V.Vector Double ]
randomColumnVectors n = Prelude.map (\vs -> M.fromLists (bundle 1 vs )) (bundle n vVals)


checkLinearSystemSolution :: M.Matrix V.Vector Double
                          -> M.Matrix V.Vector Double -> Double
                          -- -> (M.Matrix V.Vector Double, M.Matrix V.Vector Double)
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
    else error "unequal matrix dimensions is maxAbsDiff"


-- Return the largest difference (by relative value) between
-- two matrices.
--
maxRelDiff :: M.Matrix V.Vector Double
           -> M.Matrix V.Vector Double
           -> Double
maxRelDiff (Matrix m n _ _ v) (Matrix m' n' _ _ v') =
    if m == m' && n == n'
    then V.maximum $ V.zipWith (\y y' -> abs (2 * (y - y') / (abs y + abs y'))) v v'
    else error "unequal matrix dimensions is maxRelDiff"


matrixSpecRelDiff :: Double
matrixSpecRelDiff = 1.0e-9

luCheck10 :: M.Matrix V.Vector Double
luCheck10 = fromLists [[1.0], [2.0], [3.0]]

luCheck10' :: M.Matrix V.Vector Double
luCheck10' = fromLists [[1.0], [2.0], [3.0]]


-- The tests themselves
--
spec :: Spec
spec = do
    describe "Numeric.LinearAlgebra.LUSolve.luSolve" $ do
        it "random 10 * 10 matrix is solved" $ do
            maxRelDiff luCheck10 luCheck10'
                `shouldSatisfy` (< matrixSpecRelDiff)



main :: IO ()
main = hspec spec
