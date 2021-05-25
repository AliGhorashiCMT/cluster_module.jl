using Pkg;
ENV["PYTHON"] = ""
Pkg.build("PyCall")
using Test, PyCall, cluster_module

@testset "Cluster_Module" begin

end

@testset "Rotating Ions" begin
    #Check rotations about z
    @test rotateinz([1, 0, 0], π) ≈ [-1, 0, 0]
    @test rotateion([1, 0, 0], [0, 0, 0], [0, 1, 0], π) ≈ [-1, 0, 0]
    @test findperp([0, 0, 0], [1, 0, 0], [0, 1, 0]) ≈ [0, 0, 1]
    a, b, c = rand(3), rand(3), rand(3)
    @test isapprox(sum(findperp(a, b, c).*(a-b)), 0, atol=1e-3)
end

@testset "Obtaining angles and distances" begin
    @test getangle([0, 0, 0], [1, 0, 0], [-1/2, sqrt(3)/2, 0])*180/π ≈ 120    
    @test getdist([0, 0, 0], [√3/2, 1/2, 0]) ≈ 1
end

@testset "Parsing JDFTX input files" begin
    @test ionparser("ion C 0 0.3333333 0 1") ≈ [0, 0.3333333, 0]
end