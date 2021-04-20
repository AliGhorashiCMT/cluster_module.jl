#================== Methods to create cluster types from JDFTX formatted lattice and ionic positions ==================#
function ionlattice2central(ionfile::String, latticefile::String; superlattice::Vector{<:Integer}=[0, 0, 0])
    ionic_positions = np.loadtxt(ionfile, usecols=[2, 3, 4])
    initialnumions = np.shape(ionic_positions)[1]
    lattice_vectors = np.loadtxt(latticefile, skiprows=1, usecols=[0, 1, 2])
    ##Deal with superlattice 
    superx, supery, superz = superlattice
    for x in 0:superx
        for y in 0:supery
            for z in 0:superz
                for i in 1:initialnumions
                    ionic_positions = [ionic_positions; reshape(ionic_positions[i, :] + [x, y, z], (1, 3))]
                end
            end
        end
    end
    finalnumions = np.shape(ionic_positions)[1]
    iontypes = Vector{String}()
    for ionid in readlines(ionfile)
        push!(iontypes, string(split(ionid)[2]))
    end
    println(iontypes)
    cartesian_ioncoords = Vector{Vector{Float64}}()
    for ion in 1:finalnumions
        push!(cartesian_ioncoords, round.(lattice_vectors*ionic_positions[ion, :], digits=3))
    end
    #println(cartesian_ioncoords)
    println(finalnumions)

    appropriateformationpos = Vector{Vector{Any}}()
    for i in initialnumions+1:finalnumions
        try
            push!(appropriateformationpos, ["ion", iontypes[i%initialnumions], cartesian_ioncoords[i]..., 1])
        catch BoundsError
            push!(appropriateformationpos, ["ion", iontypes[9], cartesian_ioncoords[i]..., 1])
        end
    end
    #println(length(appropriateformationpos))
    [println(appropriate) for appropriate in appropriateformationpos]
    make_xsf(appropriateformationpos, lattice = [40 0 0 "\\"; 0 40 0 "\\"; 0 0 40 "\\"])
    #TODO make sure that ending backslashes don't break this coord_new_center
end