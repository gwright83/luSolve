{-# OPTIONS_GHC -fno-warn-type-defaults #-}
--
-- test/Haddock.hs
--
-- Check coverage of haddock documentation.
--
module Main (main) where

import           Data.List      (genericLength)
import           Data.Maybe     (catMaybes)
import           System.Exit    (exitFailure, exitSuccess)
import           System.Process (readProcessWithExitCode)
import           Text.Regex     (matchRegex, mkRegex)

average :: (Fractional a, Real b) => [b] -> a
average xs = realToFrac (sum xs) / genericLength xs

expected :: Fractional a => a
expected = 90

-- Note that stack sends its informational output to stderr, not
-- stdout.  Go figure.
--
main :: IO ()
main = do
    (_, _, stderrOutput)  <- readProcessWithExitCode "stack" ["haddock"] ""
    if average (match stderrOutput) >= expected
        then exitSuccess
        else putStr stderrOutput >> exitFailure

match :: String -> [Int]
match = fmap read . concat . catMaybes . fmap (matchRegex pattern) . lines
  where
    pattern = mkRegex "^ *([0-9]*)% "
