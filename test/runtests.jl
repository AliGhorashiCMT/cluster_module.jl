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

@testset "Obtaining angles and distances" begin
    getangle([0, 0, 0], [1, 0, 0], [-1/2, sqrt(3)/2, 0])*180/π ≈ 120    
    getdist([0, 0, 0], [√3/2, 1/2, 0]) ≈ 1
end

@test "Parsing JDFTX input files" begin
    ionparser("ion C 0 0.3333333 0 1") ≈ [0, 0.3333333, 0]
end