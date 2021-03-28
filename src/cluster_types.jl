abstract type cluster end 

mutable struct centralcluster<:cluster
    central_atom::Tuple{<:AbstractString, <:Vector{<:Real}}
    branch_atoms::Union{Nothing, Vector{<:Tuple{<:AbstractString, <:Vector{<:Real}}}}
end

mutable struct branchcluster<:cluster
    split_off_atom::Tuple{<:AbstractString, <:Vector{<:Real}}
    branch_atoms::Union{Nothing, Vector{<:Tuple{<:AbstractString, <:Vector{<:Real}}}} #We can also make a branch that is just one atom
end

centralcluster(central_atom::Tuple{<:AbstractString, <:Vector{<:Real}}) = centralcluster(central_atom, nothing)
branchcluster(central_atom::Tuple{<:AbstractString, <:Vector{<:Real}}) = branchcluster(central_atom, nothing)

"String macro to easily make central clusters"
macro cc_str(p)
    stringed_p = string.(split(p))
    if length(stringed_p) < 5
        error("Incorrect usage of cc string macro")
    end
    return make_central_type(stringed_p[1:2]..., parse(Int, stringed_p[3]), parse.(Float64, stringed_p[4:5])...)
end

"String macro to easily make branch clusters"
macro bc_str(p)
    stringed_p = string.(split(p))
    if length(stringed_p) < 5
        error("Incorrect usage of bc string macro")
    end
    translation = 0
    θ₁ = 0 
    θ₂ = 0
    try
        translation = parse(Float64, stringed_p[6])
        θ₁ = parse(Float64, stringed_p[7])
        θ₂ = parse(Float64, stringed_p[8])
    catch(BoundsError)

    end
    return make_branch_type(stringed_p[1:2]..., parse(Int, stringed_p[3]), parse.(Float64, stringed_p[4:5])..., translation=translation, θ₁=θ₁, θ₂=θ₂)
end

function make_central_type(middle_atom::AbstractString, linked_atom::AbstractString, num_links::Int, θ::Real, d::Real, Radian::Bool=false)

    if Radian == false ##Check if angles given in radian or not 
        θ *=  π/180
    end

    first_pos = [middle_atom]  
    second_pos = [linked_atom]  

    linked_pos = Vector{Tuple{String, Vector{Real}}}()

    θ2 = 2π/(num_links-1) # Position atoms symmetrically about the central atom
    rotation_mat = [cos(θ2) sin(θ2) 0 ; -sin(θ2) cos(θ2) 0; 0 0 1] #Rotation about z axis conventionally
    generate_pos =  [d*sin(π-θ), 0, -d*cos(π-θ)] #First position of the linked atoms

    push!(linked_pos, (linked_atom, [0, 0, d]))
    for n in 0:num_links-2
        push!(linked_pos, (linked_atom, [(rotation_mat^n)*generate_pos...]))
    end
    return centralcluster((middle_atom, [0, 0, 0]), linked_pos)

end

function make_branch_type(middle_atom::AbstractString, linked_atom::AbstractString, num_links::Int, θ::Real, d::Real; translation::Real=3, θ₁::Real=0, θ₂::Real=0, Radian::Bool=false)

    if Radian == false ##Check if angles given in degrees or radians, and convert accordingly 
        θ₁ *= π/180
        θ₂ *= π/180
        θ *= π/180
    end
    #=The translation passed to mak_branch_type is the relative translation of the central atom of the branch with respect to 
    the branching point of the central cluster
    =#
    additional_rotation1 = [1 0 0; 0 cos(θ₁) -sin(θ₁); 0 sin(θ₁) cos(θ₁)  ] ##Rotation about x axis
    additional_rotation2 = [cos(θ₂) 0 sin(θ₂); 0 1 0; -sin(θ₂) 0 cos(θ₂)  ] ##Rotation about z axis

    first_pos = [(additional_rotation2*additional_rotation1*[0, 0, -translation])...]
    linked_pos = Vector{Tuple{String, Vector{Real}}}()
    rotation_mat = [cos(2π/num_links) sin(2π/num_links) 0 ; - sin(2π/num_links) cos(2π/num_links) 0; 0 0 1] #Rotate the initial generating position about the z axis
    generate_pos =  [d*sin(π-θ), 0, -d*cos(π-θ)-translation] #The generating position

    for n in 0:num_links-1
        push!(linked_pos, (linked_atom, [additional_rotation2*additional_rotation1*rotation_mat^n*generate_pos...]))
    end
    return branchcluster((middle_atom, first_pos), linked_pos)
end


function move_branch(central_cluster::centralcluster, branch_cluster::branchcluster, which_atom::Int)

    ##Returns the moved branch with respect to the central cluster

    coord_new_center = central_cluster.branch_atoms[which_atom-1][2];
    
    d = sqrt(sum(coord_new_center.^2));
    α, θ = find_α_θ(coord_new_center);

    rotation_mat = [ cos(θ) -sin(θ)cos(α) sin(α)sin(θ); sin(θ) cos(θ)cos(α) -sin(α)cos(θ); 0 sin(α) cos(α)];
    moved_branch_ids = Vector{Tuple{AbstractString, Vector{Real}}}();

    for atom_id in branch_cluster.branch_atoms
        push!(moved_branch_ids,  (atom_id[1], [(rotation_mat*(atom_id[2]-[0, 0, d])...)]  ) )
    end

    new_split_off_atom = (branch_cluster.split_off_atom[1], rotation_mat*(branch_cluster.split_off_atom[2]-[0, 0, d]))
    return branchcluster(new_split_off_atom, moved_branch_ids)

end

#= 
Methods to combine clusters together 
=#

function Base.append!(branch_cluster::branchcluster, branch_cluster2::branchcluster)
    return branch_cluster
end

function Base.append!(central_cluster::centralcluster, branch_cluster::branchcluster)
    temp_branch = branchcluster(branch_cluster.split_off_atom, copy(branch_cluster.branch_atoms))
    temp_central = centralcluster(central_cluster.central_atom, copy(central_cluster.branch_atoms))
    #return centralcluster(temp_central.central_atom, append!(temp_branch.branch_atoms, append!(temp_central.branch_atoms, [temp_branch.split_off_atom])))
    return centralcluster(temp_central.central_atom, append!(temp_central.branch_atoms, append!(temp_branch.branch_atoms, [temp_branch.split_off_atom])))
end

function Base.append!(branch_cluster::branchcluster, central_cluster::centralcluster)
    temp_branch = branchcluster(branch_cluster.split_off_atom, copy(branch_cluster.branch_atoms))
    temp_central = centralcluster(central_cluster.central_atom, copy(central_cluster.branch_atoms))
    return centralcluster(temp_central.central_atom, append!(temp_branch.branch_atoms, append!(temp_central.branch_atoms, [temp_branch.split_off_atom])))
end

function Base.append!(a::Nothing, b::Any) #To prevent errors when appending central clusters with no branch atoms
    return b
end

function Base.append!(b::Any, a::Nothing ) #To prevent errors when appending central clusters with no branch atoms
    return b
end

function Base.copy(a::Nothing)
    return nothing
end

function +(central_cluster::centralcluster, branch_cluster::branchcluster)
    return append!(central_cluster, branch_cluster)
end

function +(branch_cluster::branchcluster, central_cluster::centralcluster)
    return append!(central_cluster, branch_cluster)
end

function Base.length(branch_cluster::branchcluster)
    return length(branch_cluster.branch_atoms) + 1
end

function Base.length(central_cluster::centralcluster)
    return length(central_cluster.branch_atoms) + 1
end

function Base.length(a::Nothing) ## To prevent errors in other length methods
    return 0 
end

function Base.show(io::IO, central_cluster::centralcluster )
    println(io, "Central atom is: $(central_cluster.central_atom[1])  at coordinate $(round.(central_cluster.central_atom[2], digits=3))" )
    if isa(central_cluster.branch_atoms, Nothing) == false
        for atom in central_cluster.branch_atoms
            println(io, "Branch atom: $(atom[1])  ", "at coordinate   $(round.(atom[2], digits=3))")
        end
    end
    println("\n\n Total Atoms: $(length(central_cluster))")
end

function Base.show(io::IO, branch_cluster::branchcluster)
    println(io, "Split off atom is: $(branch_cluster.split_off_atom[1])  at coordinate $(round.(branch_cluster.split_off_atom[2], digits=3))" )
    if isa(branch_cluster.branch_atoms, Nothing) == false
        for atom in branch_cluster.branch_atoms
            println(io, "Branch atom: $(atom[1])  ", "at coordinate   $(round.(atom[2], digits=3))")
        end
    end
    println("\n\n Total atoms: $(length(branch_cluster))")

end

function make_new_central(old_central::centralcluster, branch_cluster::branchcluster, where_attach::Vector{Int})
    for (index, atom) in enumerate(where_attach)
        old_central = append!(old_central, move_branch(old_central, branch_cluster, atom));
    end
    return old_central;
end

function make_new_central(old_central::centralcluster, branches::Vector{branchcluster})

    for (index, branch) in enumerate(branches)
        old_central = append!(old_central, branch )
    end
    return old_central 
end


function make_xsf(cluster::centralcluster; lattice::Array{<:Any, 2} = [40 0 0 "\\"; 0 40 0 "\\"; 0 0 40 "\\"])
    ##Convert to expected JDFTX ionpos format
    ionpos = Vector{Vector{Any}}()
    push!(ionpos, ["ion", cluster.central_atom[1], cluster.central_atom[2]..., 1 ])
    for i in 1:length(cluster)-1
        push!(ionpos, ["ion", cluster.branch_atoms[i][1], cluster.branch_atoms[i][2]..., 1 ])
    end
    make_xsf(ionpos, lattice)
end
