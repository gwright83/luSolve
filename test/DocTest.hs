--
-- test/DocTest.hs
--
-- Run doctest on any code examples in haddock comments.
--
module Main (main) where

import           System.FilePath.Glob (glob)
import           Test.DocTest         (doctest)

main :: IO ()
main = glob "library/**/*.hs" >>= doctest
