# cluster_module.jl

[![Build status][ci-status-img]][ci-status-url][![Coverage][codecov-img]][codecov-url]

A Julia module for creating large molecular clusters for DFT calculations by either Abinit, JDFTX or Quantum Espresso 

For plotting functionalities, JDFTX, JDFTX post-processing bash scripts, and Vesta must be installed and in the $PATH

## Installation 
```julia
pkg> add https://github.com/AliGhorashiCMT/cluster_module.jl
```
which will allow access via
```julia
julia> using cluster_module
```

[ci-status-img]:   https://github.com/AliGhorashiCMT/cluster_module.jl/workflows/CI/badge.svg
[ci-status-url]:   https://github.com/AliGhorashiCMT/cluster_module.jl/actions
[codecov-img]: https://codecov.io/gh/AliGhorashiCMT/Jcluster_module.jl/branch/main/graph/badge.svg
[codecov-url]: https://app.codecov.io/gh/AliGhorashiCMT/cluster_module.jl
