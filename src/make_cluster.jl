
abstract type cluster end 

struct centralcluster<:cluster
    central_atom::Array{<:Any, 1}
    branch_atoms::Array{<:Array{<:Any, 1}}
end

struct branchcluster<:cluster
    atoms::Array{<:Array{<:Any, 1}}
end

function make_xsf(ionpos::Array{<:Array{<:Any, 1}, 1}; lattice::Array{<:Any, 2} = [40 0 0 "\\"; 0 40 0 "\\"; 0 0 40 "\\"])
    print("making xsf")
    open("temp.in", "w") do io
        write(io, "ion-species GBRV/\$ID_pbesol.uspp")
        writedlm(io, "  ")
        write(io, "core-overlap-check none")

        writedlm(io, "  ")
        write(io, "coords-type Cartesian")
        writedlm(io, "  ")
        writedlm(io, ionpos, " ")
        write(io, "lattice\\")
        writedlm(io, " ")
        writedlm(io, lattice, " " )
    end;
    run(pipeline(`jdftx -ni temp.in `, `tee temp.out`));
    run(`createXSF temp.out temp.xsf`);
    run(`open temp.xsf`)
end


function make_xsf(ionpos::Any, lattice::Array{<:Any, 2} = [20 0 0 "\\"; 0 20 0 "\\"; 0 0 20 "\\"])
    open("temp.in", "w") do io
        write(io, "ion-species GBRV/\$ID_pbesol.uspp")
        writedlm(io, "  ")
        write(io, "core-overlap-check none")
        writedlm(io, "  ")
        write(io, "coords-type Cartesian")
        writedlm(io, "  ")
        writedlm(io, ionpos, " ")
        write(io, "lattice\\")
        writedlm(io, " ")
        writedlm(io, lattice, " " )
    end;
    run(pipeline(`jdftx -ni temp.in `, `tee temp.out`));
    run(`createXSF temp.out temp.xsf`);
    run(`open temp.xsf`)
end


function make_test_xsf()
    path_to_test_si = joinpath(@__DIR__, "../examples/si.out");
    if isfile("si-copy.out")
        rm("si-copy.out")
    end
    cp( path_to_test_si,  "si-copy.out")
    run(`createXSF si-copy.out si.xsf `)
    run(`open si.xsf`)
    rm("si-copy.out")
end

function find_α_θ(rotated_vector::Vector{<:Any})
    normalized_vector = rotated_vector/sqrt(sum(rotated_vector.^2))
    α = acos(-normalized_vector[3])
    θ = asin(-normalized_vector[1]/sin(α))
    if (sin(α)*cos(θ) ≈ normalized_vector[2]) == false
        α = -α; θ = -θ
    end
    return α, θ
end

function attach_cluster(central_cluster::Any, branch_clusters::Array{<:Any}, split_points::Vector{Int})
    complete_cluster = Any[]
    push!(complete_cluster, central_cluster...)
    for i in 1:length(branch_clusters)
        push!(complete_cluster, move_cluster(central_cluster, branch_clusters[i], split_points[i] )...)
    end
    return complete_cluster
end

function move_cluster(central_cluster::Array{<:Any}, new_cluster::Array{<:Any}, which_atom::Int)

    coord_new_center = central_cluster[which_atom][3:5]
    d = sqrt(sum(coord_new_center.^2))
    α, θ = find_α_θ(coord_new_center)

    rotation_mat = [ cos(θ) -sin(θ)cos(α) sin(α)sin(θ); sin(θ) cos(θ)cos(α) -sin(α)cos(θ); 0 sin(α) cos(α)]
    new_cluster1 = []

    for atom in new_cluster
        push!(new_cluster1, [atom[1], atom[2], (rotation_mat*(atom[3:5]-[0, 0, d]))..., 1] )
    end
    return new_cluster1
end

function make_arbitrary_cluster(middle_atom::AbstractString, linked_atom::AbstractString, num_links::Int, θ::Real, d::Real; translation::Real=3, θₓ::Real=0, θ₂::Real=0)

    additional_rotationx = [1 0 0; 0 cos(θₓ) -sin(θₓ); 0 sin(θₓ) cos(θₓ)  ]
    additional_rotation2 = [cos(θ₂) 0 sin(θ₂); 0 1 0; -sin(θ₂) 0 cos(θ₂)  ]

    first_pos = ["ion", middle_atom, (additional_rotation2*additional_rotationx*[0, 0, -translation])..., 1]  
    linked_pos = Array{Array{Any, 1}, 1}(undef, num_links+1)
    rotation_mat = [cos(2π/num_links) sin(2π/num_links) 0 ; - sin(2π/num_links) cos(2π/num_links) 0; 0 0 1]
    generate_pos =  [d*sin(π-θ), 0, -d*cos(π-θ)-translation]
    linked_pos[1] = first_pos
    for n in 0:num_links-1
        linked_pos[n+2] =["ion", linked_atom, (additional_rotation2*additional_rotationx*((rotation_mat^n)*generate_pos))..., 1]
    end
    return linked_pos

end

function make_central_cluster(middle_atom::AbstractString, linked_atom::AbstractString, num_links::Int, θ::Real, d::Real)

    first_pos = ["ion", middle_atom, 0, 0, 0, 1]  
    second_pos = ["ion", linked_atom, 0, 0, d, 1]  

    linked_pos = Array{Array{Any, 1}, 1}(undef, num_links+1)
    θ2 = 2π/(num_links-1)
    rotation_mat = [cos(θ2) sin(θ2) 0 ; -sin(θ2) cos(θ2) 0; 0 0 1]
    generate_pos =  [d*sin(π-θ), 0, -d*cos(π-θ)]

    linked_pos[1] = first_pos
    linked_pos[2] = second_pos
    for n in 0:num_links-2
        linked_pos[n+3] = ["ion", linked_atom, (rotation_mat^n)*generate_pos..., 1]
    end
    return linked_pos

end
function make_alpha_cluster()
    lat_vectors = [4.604*1.889 -2.3020000457999998*1.889 0 "\\"; 0 3.987*1.889 0 "\\"; 0 0 5.207*1.889 "\\"  ]
    ionpos =   [["ion", "Si", 0.4436617824484789 ,-0.0000000000000000 , 0.3333333429999996, 1],
    ["ion", "Si", -0.0000000000000000,  0.4436617824484789 , 0.6666666870000029, 1],
     ["ion", "Si", 0.5563382175515210,  0.5563382175515210 ,-0.0000000000000000, 1],
     ["ion", "O", 0.3926661416221499 , 0.3062177364999842  ,0.2428214976299141, 1],
     ["ion", "O", 0.6937822635000156 , 0.0864484051221655  ,0.5761548406299137, 1],
     ["ion", "O", 0.9135515948778347 , 0.6073338583778505  ,0.9094881546299145, 1],
     ["ion", "O", 0.3062177364999842 , 0.3926661416221499  ,0.7571785323700884, 1],
     ["ion", "O", 0.0864484051221655 , 0.6937822635000156  ,0.4238451593700863, 1],
     ["ion", "O", 0.6073338583778505  ,0.9135515948778347  ,0.0905118383700884, 1]]
     return ionpos, lat_vectors
end