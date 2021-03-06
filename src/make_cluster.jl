
function make_xsf(ionpos::Array{<:Array{<:Any, 1}, 1}; lattice::Array{<:Any, 2} = [1000 0 0 "\\"; 0 1000 0 "\\"; 0 0 1000 "\\"])
    print("making xsf")
    open("temp.in", "w") do io
        write(io, "ion-species GBRV/\$ID_pbesol.uspp")
        writedlm(io, "  ")
        #write(io, "core-overlap-check none")
        write(io, "  ")
        #write(io, "\ncoulomb-interaction Isolated")
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

function run_jdftx(;filename::AbstractString="temp.in", parallel::Bool=false)
    if parallel 
        run(pipeline(`mpirun -n 4 jdftx -i $(filename)`,  `tee $(filename[1:end-3]*".out")` ))  
    else 
        run(pipeline(`jdftx -i $(filename)`,  `tee $(filename[1:end-3]*".out")` ))  
    end
end

function run_jdftx_ni(;filename::AbstractString="temp.in", parallel::Bool=false)
    if parallel 
        run(pipeline(`mpirun -n 4 jdftx -ni $(filename)`,  `tee $(filename[1:end-3]*".out")` ))  
    else
        run(pipeline(`jdftx -ni $(filename)`,  `tee $(filename[1:end-3]*".out")` ))  
    end
end

function make_xsf(ionpos::Any, lattice::Array{<:Any, 2} = [1000 0 0 "\\"; 0 1000 0 "\\"; 0 0 1000 "\\"])
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

function find_??_??(rotated_vector::Vector{<:Any})
    normalized_vector = rotated_vector/sqrt(sum(rotated_vector.^2))
    ?? = acos(-normalized_vector[3])
    ?? = asin(-normalized_vector[1]/sin(??))
    if (sin(??)*cos(??) ??? normalized_vector[2]) == false
        ?? = -??; ?? = -??
    end
    if isequal(??, NaN)
        ?? = 0 
    end

    if isequal(??, NaN)
        ?? = 0 
    end

    rotation_mat = [ cos(??) -sin(??)cos(??) sin(??)sin(??); sin(??) cos(??)cos(??) -sin(??)cos(??); 0 sin(??) cos(??)]

    if (rotation_mat*[0, 0, -1] ??? normalized_vector) == false ##Throws an error if for some reason the rotation did something unexpected
        error("find_??_?? didn't move cluster correctly")
    end
    return ??, ??

end

function attach_cluster(central_cluster::Any, branch_clusters::Array{<:Any}, split_points::Vector{Int})
    complete_cluster = Any[]
    push!(complete_cluster, central_cluster...)
    for i in 1:length(branch_clusters)
        push!(complete_cluster, move_cluster(central_cluster, branch_clusters[i], split_points[i] )...)
    end
    return complete_cluster
end

function attach_clusters(central_cluster::Any, branch_cluster::Array{Array{Any,1},1}, split_points::Vector{Int})
    return attach_cluster(central_cluster, repeat(branch_cluster, length(split_points)), split_points)
end

function move_cluster(central_cluster::Array{<:Any}, new_cluster::Array{<:Any}, which_atom::Int)

    coord_new_center = central_cluster[which_atom][3:5]
    d = sqrt(sum(coord_new_center.^2))
    ??, ?? = find_??_??(coord_new_center)

    rotation_mat = [ cos(??) -sin(??)cos(??) sin(??)sin(??); sin(??) cos(??)cos(??) -sin(??)cos(??); 0 sin(??) cos(??)]
    new_cluster1 = []

    for atom in new_cluster
        push!(new_cluster1, [atom[1], atom[2], (rotation_mat*(atom[3:5]-[0, 0, d]))..., 1] )
    end
    return new_cluster1
end

function make_arbitrary_cluster(middle_atom::AbstractString, linked_atom::AbstractString, num_links::Int, ??::Real, d::Real; translation::Real=3, ?????::Real=0, ?????::Real=0)

    additional_rotationx = [1 0 0; 0 cos(?????) -sin(?????); 0 sin(?????) cos(?????)  ]
    additional_rotation2 = [cos(?????) 0 sin(?????); 0 1 0; -sin(?????) 0 cos(?????)  ]

    first_pos = ["ion", middle_atom, (additional_rotation2*additional_rotationx*[0, 0, -translation])..., 1]  
    linked_pos = Array{Array{Any, 1}, 1}(undef, num_links+1)
    rotation_mat = [cos(2??/num_links) sin(2??/num_links) 0 ; - sin(2??/num_links) cos(2??/num_links) 0; 0 0 1]
    generate_pos =  [d*sin(??-??), 0, -d*cos(??-??)-translation]
    linked_pos[1] = first_pos
    for n in 0:num_links-1
        linked_pos[n+2] =["ion", linked_atom, (additional_rotation2*additional_rotationx*((rotation_mat^n)*generate_pos))..., 1]
    end
    return linked_pos

end

function make_central_cluster(middle_atom::AbstractString, linked_atom::AbstractString, num_links::Int, ??::Real, d::Real)

    first_pos = ["ion", middle_atom, 0, 0, 0, 1]  
    second_pos = ["ion", linked_atom, 0, 0, d, 1]  

    linked_pos = Array{Array{Any, 1}, 1}(undef, num_links+1)
    ??2 = 2??/(num_links-1)
    rotation_mat = [cos(??2) sin(??2) 0 ; -sin(??2) cos(??2) 0; 0 0 1]
    generate_pos =  [d*sin(??-??), 0, -d*cos(??-??)]

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