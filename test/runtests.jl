using Pkg;
ENV["PYTHON"] = ""
Pkg.build("PyCall")
using Test, PyCall, cluster_module

@testset "Cluster_Module" begin

end

@testset "Rotating Ions" begin
    #Check rotations about z
    rotateinz([1, 0, 0], π) ≈ [-1, 0, 0]
    rotateion([1, 0, 0], [0, 0, 0], [0, 1, 0], π) ≈ [-1, 0, 0]
    findperp([0, 0, 0], [1, 0, 0], [0, 1, 0]) ≈ [0, 0, 1]
    a, b, c = rand(3), rand(3), rand(3)
    sum(findperp(a, b, c).*(a-b))
end