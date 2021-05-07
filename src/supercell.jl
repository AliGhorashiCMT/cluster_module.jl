"""
$(TYPEDSIGNATURES)
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


