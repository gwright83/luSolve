--
-- Benchmarks for LUSolve.
--

module LUSolveBenchmark (
    benchmarks
    ) where

import           Criterion                     (Benchmark, bench, nf)

import           Numeric.LinearAlgebra.LUSolve (luFactor)

import qualified Data.Matrix.Dense.Generic     as M
import qualified Data.Vector.Unboxed           as V
import           System.Random

bundle :: Int -> [ a ] -> [[ a ]]
bundle _ [] = []
bundle n xs = take n xs : bundle n (drop n xs)

mVals :: [ Double ]
mVals = randoms (mkStdGen 1)


randomSquareMatrices :: Int -> [ M.Matrix V.Vector Double ]
randomSquareMatrices n = Prelude.map (\vs -> M.fromLists (bundle n vs)) (bundle (n * n) mVals)

_vVals :: [ Double ]
_vVals = randoms (mkStdGen 2)

_randomColumnVectors :: Int -> [ M.Matrix V.Vector Double ]
_randomColumnVectors n = Prelude.map (\vs -> M.fromLists (bundle 1 vs )) (bundle n _vVals)

runLUFactor :: Int -> M.Matrix V.Vector Double
runLUFactor n = (\(x, _, _) -> x) $ luFactor $ head $ randomSquareMatrices n


benchmarks :: [Benchmark]
benchmarks =
    [ bench "luFactor 100 x 100 matrix"   $ nf runLUFactor 100
    , bench "luFactor 500 x 500 matrix"   $ nf runLUFactor 500
    , bench "luFactor 1000 x 1000 matrix" $ nf runLUFactor 1000
    ]
