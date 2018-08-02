luSolve
-------

A pure Haskell implementation of LU decomposition, along with a solver
for linear systems in the form Ax = b.  Partial pivoting is used to
avoid numerical failure on otherwise well-conditioned input matrices.

To see something happen, type
```
    $ stack test
```
