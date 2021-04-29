using Pkg;
ENV["PYTHON"] = ""
Pkg.build("PyCall")
using Test, PyCall, cluster_module

@testset "Cluster_Module" begin

end
