--
-- Benchmarks for LUSolve.
--

module LUSolveBenchmark (
    benchmarks
    ) where

import           Criterion                     (Benchmark, bench, bgroup, env,
                                                nf)

import           Numeric.LinearAlgebra.LUSolve (luFactor)

import qualified Data.Matrix.Dense.Generic     as M
import qualified Data.Vector.Unboxed           as V
import           System.Random

bundle :: Int -> [ a ] -> [[ a ]]
bundle _ [] = []
bundle n xs = take n xs : bundle n (drop n xs)

type Mat = M.Matrix V.Vector Double

runLUFactor :: Mat -> Mat
runLUFactor = (\(x, _, _) -> x) . luFactor

setupEnv :: IO (Mat, Mat, Mat)
setupEnv = do
  let mVals = randoms (mkStdGen 1) -- not a top level CAF, so will be GC'd promptly
      randomSquareMatrices :: Int -> [ Mat ]
      randomSquareMatrices n = Prelude.map (\vs -> M.fromLists (bundle n vs)) (bundle (n * n) mVals)
      m100:_ = randomSquareMatrices 100
      m500:_ = randomSquareMatrices 500
      m1000:_= randomSquareMatrices 1000
  return (m100, m500, m1000)

benchmarks :: Benchmark
benchmarks = env setupEnv $ \ ~(m100, m500, m1000) ->
  bgroup "luFactor"
    [ bench "100 x 100"   $ nf runLUFactor m100
    , bench "500 x 500"   $ nf runLUFactor m500
    , bench "1000 x 1000" $ nf runLUFactor m1000
    ]
