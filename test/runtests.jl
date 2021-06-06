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

@testset "Creating Supercells" begin
    #We use a model for alpha quartz. Note that the ionpos files are interpreted as being in the lattice basis
    #Lattice is also written in Bohr
    #We make a supercell of mult 2 2 2
    i, l = createsupercell("alpha.ionpos", "alpha.lattice", [2, 2, 2], writefile="alpha.supercell.in")
    @test length(i) == 9*8 #Check that correct number of ions are being produced
end

@testset "Creating a cluster from a supercell" begin
    #We make a cluster model for a self trapped hole defect in alpha quartz silica
    i, l = createsupercellcluster("alpha.ionpos", "alpha.lattice", [3, 3, 2], writefile="alpha.cluster.in", removecriteria=(157, 10.2))
    @test length(i) < 3*3*2*9
    i, l = createsupercellcluster("alpha.ionpos", "alpha.lattice", [3, 3, 2], writefile="alpha.cluster.in")
    @test length(i) == 162
end

@testset "Extending Bonds Tests" begin
    newion = extendbond([1/4, 1/4, 1/4], [3/4, 3/4, 3/4], 5)
    isapprox(sqrt(sum((newion-[0.25, 0.25, 0.25]).^2 )), 5, atol=1e-3)
end

