--
-- Benchmark driver for luSolve.
--

module Main (main) where

import           Criterion.Main   (defaultMain)
import qualified LUSolveBenchmark


main :: IO ()
main = defaultMain [ LUSolveBenchmark.benchmarks ]
