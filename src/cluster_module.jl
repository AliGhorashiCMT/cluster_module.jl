module cluster_module

using DelimitedFiles, PyCall, DocStringExtensions, LinearAlgebra

import Base: + 

const np = PyNULL()
const ase = PyNULL()
const ase_atoms = PyNULL()
const atoms = PyNULL()
export np
export ase
export ase_atoms
export atoms
# ---------------------------------------------------------------------------------------- #
function __init__()
    copy!(np, pyimport_conda("numpy", "numpy"))
    copy!(ase, pyimport_conda("ase", "ase", "conda-forge"))
    copy!(ase_atoms, pyimport_conda("ase.atoms", "ase", "conda-forge"))
    copy!(atoms, ase_atoms.Atoms )
end

include("make_cluster.jl")
export make_xsf
export make_test_xsf, find_α_θ
export make_alpha_cluster, make_central_cluster, attach_cluster, make_arbitrary_cluster, run_jdftx, run_jdftx_ni
export cluster, centralcluster, branchcluster

include("cluster_types.jl")
export make_central_type, make_branch_type, move_branch, make_new_central, @cc_str, @bc_str

include("clusterfromionpos.jl")
export ionlattice2central, deleteions

include("rotateions.jl")
export rotateion, rotateions

include("getangles.jl")
export getangle, tocartesian, tocartesians, getdist, extendbond

include("supercell.jl")
export createsupercell, createsupercellcluster, replaceions
end # module
