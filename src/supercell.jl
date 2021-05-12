"""
$(TYPEDSIGNATURES)

`ionpos` : Name of JDFTX ionic positions file 

`lattice` : Name of JDFTX lattice file 

`supermults` : Vector denoting the supercell multiplicity

`writefile` : If nothing indicates that the supercell is not to be written to a file. If a string, is taken to be the name of a file to which the supercell will be written.
"""
function createsupercell(ionpos::String, lattice::String,  supermults::Vector{<:Real}; writefile::Union{Nothing, String}=nothing)
    m1, m2, m3 = supermults
    readions = Vector{Vector{Real}}()
    ionids = Vector{Vector{Any}}()
    supercellions = Vector{Vector{Any}}()
    latt = Array{Float64, 2}(undef, (3, 3))
    counter = 1
    for line in readlines(lattice)
        try
            latt[counter, :] = parse.(Float64, String.(split(replace(line, "\\"=>""))[1:3])) #Parses the lattice file
            counter += 1
        catch
        end
    end
    for line in readlines(ionpos)
        push!(readions, parse.(Float64, String.(split(line)[3:5])))
        push!(ionids, String.(split(line)[1:2]))
    end
    for (ionid, ion) in zip(ionids, readions)
        for (n1, n2, n3) in Tuple.(CartesianIndices(rand(m1, m2, m3)))
            k1, k2, k3 = (n1-1)/m1, (n2-1)/m2, (n3-1)/m3
            push!(supercellions, [ionid..., (ion./[m1, m2, m3]+[k1, k2, k3])..., 1 ])
        end
    end
    supercelllattice = latt .* [m1 m2 m3]
    if !isnothing(writefile)
        open(writefile, "w") do io 
            write(io, "ion-species GBRV/\$ID_pbesol.uspp \n"  )
            write(io, "coords-type Lattice\n")
            for ion in supercellions
                for ionparam in ion 
                    write(io, string(ionparam), "\t" )
                end
                write(io, "\n")
            end
            write(io, "\n\n")
            write(io, "lattice\\ \n")
            for (index, row) in enumerate(eachrow(supercelllattice))
                write(io, [" "*string(r) for r in collect(row) ]...)
                index < 3 ? write(io, " \\ \n") : write(io, "\n")
            end
        end
    end
    return supercellions, supercelllattice
end 


"""
$(TYPEDSIGNATURES)

`ionpos` : Name of JDFTX ionic positions file

`lattice` : Name of JDFTX lattice file

`supermults` : Vector indicating supercell multiplicity

`writefile` : If nothing- indicates cluster will not be written. If a string, indicates name of file for cluster to be written to. 

`removecriteria` : Criteria for removing atoms from the cluster. First component of tuple is the index of the atom from which distances will be measured. 
                   Second component of tuple is the maximum distance from the aforementioned atom- this is the criteria for removing atoms from the cluster. 

`removeindividual` : Indices of individual atoms to be removed.
"""
function createsupercellcluster(ionpos::String, lattice::String,  supermults::Vector{<:Real}; writefile::Union{Nothing, String}=nothing,
    removecriteria::Union{Tuple{<:Integer, <:Real}, Nothing}=nothing, removeindividual::Union{Vector{<:Integer}, Nothing}=nothing)
    m1, m2, m3 = supermults
    readions = Vector{Vector{Real}}()
    ionids = Vector{Vector{Any}}()
    supercellions = Vector{Vector{Any}}()
    supercellions2 = Vector{Vector{Any}}()
    supercellions3 = Vector{Vector{Any}}()
    latt = Array{Float64, 2}(undef, (3, 3))
    counter = 1
    for line in readlines(lattice)
        try
            latt[counter, :] = parse.(Float64, String.(split(replace(line, "\\"=>""))[1:3])) #Parses the lattice file
            counter += 1
        catch
        end
    end
    supercelllattice = latt .* [m1 m2 m3]
    supercelllatticeCART = [50 0 0; 0 50 0; 0 0 50]
    for line in readlines(ionpos)
        push!(readions, parse.(Float64, String.(split(line)[3:5])))
        push!(ionids, String.(split(line)[1:2]))
    end
    for (ionid, ion) in zip(ionids, readions)
        for (n1, n2, n3) in Tuple.(CartesianIndices(rand(m1, m2, m3)))
            k1, k2, k3 = (n1-1)/m1, (n2-1)/m2, (n3-1)/m3
            push!(supercellions2, [ionid..., (tocartesian(ion./[m1, m2, m3]+[k1, k2, k3], supercelllattice))..., 1 ])
        end
    end
    if !isnothing(removecriteria)
        baseatomvec = supercellions2[removecriteria[1]][3:5]
        for supercellion in supercellions2
            if sqrt(sum((supercellion[3:5]-baseatomvec).^2))<removecriteria[2]
                push!(supercellions3, supercellion)
            end
        end
    elseif isnothing(removecriteria)
        supercellions3 = supercellions2
    end
    if !isnothing(removeindividual)
        for (index, ion) in enumerate(supercellions3)
            index ∉ removeindividual && push!(supercellions, ion)
        end
    elseif isnothing(removeindividual)
        supercellions = supercellions3
    end
    if !isnothing(writefile)
        open(writefile, "w") do io 
            write(io, "ion-species GBRV/\$ID_pbesol.uspp \n"  )
            write(io, "coords-type Cartesian\n")
            for ion in supercellions
                for ionparam in ion 
                    write(io, string(ionparam), "\t" )
                end
                write(io, "\n")
            end
            write(io, "\n\n")
            write(io, "lattice\\ \n")
            for (index, row) in enumerate(eachrow(supercelllatticeCART))
                write(io, [" "*string(r) for r in collect(row) ]...)
                index < 3 ? write(io, " \\ \n") : write(io, "\n")
            end
        end
    end
    return supercellions, supercelllattice
end 

"""
$(TYPEDSIGNATURES)

`ionpos` : Name of the ionic positions file

`replacement` : Name of the replacement ion

`indxs` : Indices of ions to be replaced

`newfile` : A new file to which the new ion IDs will be written (if nothing- the ions are written to ionpos)
"""
function replaceions(ionpos::String, replacement::String, indxs::Vector{<:Integer}; newfile::Union{String, Nothing}=nothing)
    isnothing(newfile) && (newfile = ionpos)
    ionposes = Vector{Vector{Any}}()
    idx = 1
    for line in readlines(ionpos)
        isempty(line) && continue
        String.(split(line))[1] == "ion" || continue
        idx ∈ indxs && (push!(ionposes, ["ion", replacement, ionparser(line)..., 1]); idx+=1; continue)
        idx ∉ indxs && (push!(ionposes, [String.(split(line))[1:2]..., ionparser(line)..., 1 ]); idx+=1; continue)
    end
    open(newfile, "w") do io 
        for ion in ionposes
            write(io, [string(ionparam)*" " for ionparam in ion]..., "\n")
        end
    end
    return ionposes
end

"""
$(TYPEDSIGNATURES)

This method is provided for user ease in parsing JDFTX ionic positions files. 
The typical JDFTX ionic position is given in a form like: 

`ion C 0 1/3 2/3 0`

This method parses the line from the third to the fifth components and returns the associated vector of Floats.


`line` : A string of a JDFTX ionic position 

"""
function ionparser(line::String)
    return parse.(Float64, String.(split(line)[3:5]))
end