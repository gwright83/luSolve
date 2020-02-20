Changelog for luSolve
---------------------

# Version 0.6

On MacOS, use MacPorts llvm version 6.0 when generating code.  This gives improved
numerical performance compared to the standard GHC code generator.

If MacPorts llvm is not available, it is necessary to edit the ``luSolve.cabal``
file to set the ``useLLVM`` flag to ``false``.

## Unreleased changes
