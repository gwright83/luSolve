[![pipeline status](https://gitlab.18clay.com/software/luSolve/badges/streamlined/pipeline.svg)](https://gitlab.18clay.com/software/luSolve/commits/streamlined)

luSolve
-------

A pure Haskell implementation of LU decomposition, along with a solver
for linear systems in the form Ax = b.  Partial pivoting is used to
avoid numerical failure on otherwise well-conditioned input matrices.

The algorithm used is Gustavson's recursive LU decomposition, see
**Gustavson, F.G., "Recursion leads to automatic variable blocking for
dense linear algebra algorithms," IBM J. Res. Dev., 41, pp. 737-756 (1997).**
