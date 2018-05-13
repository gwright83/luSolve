--
-- Benchmark driver for luSolve.
--

module Main (main) where

import           Criterion.Main   (bgroup, defaultMain)
import qualified LUSolveBenchmark


main :: IO ()
main = defaultMain
    [ bgroup "LUSolve" LUSolveBenchmarks.benchmarks
    ]
