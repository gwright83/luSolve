--
-- Spec file for luSolve.hs
--
-- Tests the exposed LU decomposition and linear system solver.
--
--
module Numeric.SpecialFunction.LUSolveSpec where

import           Test.Hspec

import           Numeric.LinearAlgebra.LUSolve




-- Return the largest difference (by absolute value) between
-- two lists of Doubles.
--
maxAbsDiff :: [ Double ] -> [ Double ] -> Double
maxAbsDiff f f' = maximum $ zipWith (\y y' -> abs (y - y')) f f'


-- Return the largest difference (by relative value) between
-- two lists of Doubles.
--
maxRelDiff :: [ Double ] -> [ Double ] -> Double
maxRelDiff f f' = maximum $ zipWith (\y y' -> abs ((y - y') / y') ) f f'


-- The tests themselves
--
spec :: Spec
spec = do
    describe "Numeric.SpecFunction.BesselK0.besselK0" $ do
        it "agrees with Gnu Scientific Library gsl_sf_bessel_K0" $ do
            (besselArgs, besselK0Vals) <- readTable besselK0TableName
            maxRelDiff (map besselK0 besselArgs) besselK0Vals
                `shouldSatisfy` (< besselSpecRelDiff)



main :: IO ()
main = hspec spec
