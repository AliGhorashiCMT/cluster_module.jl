"""
$(TYPEDSIGNATURES)

Returns the angle between the vectors formed by `v1` and `v2` where 

`v1` = `branch1` - `centeratom`

and 

`v2` = `branch2` - `centeratom`
"""
function getangle(centeratom::Vector{<:Real}, branch1::Vector{<:Real}, branch2::Vector{<:Real})
    #Take dot product of the two displacement vectors and divide by the lengths 
    #Then take arccos of answer
    dist1 = branch1-centeratom
    dist2 = branch2-centeratom
    println(sqrt(sum(dist1.^2))*0.529, " ", sqrt(sum(dist2.^2))*0.529)
    println(round(acos(sum(dist1.*dist2)/(sqrt(sum(dist1.^2))*sqrt(sum(dist2.^2))))/π, digits=3), " in units of pi" )
    println(round(acos(sum(dist1.*dist2)/(sqrt(sum(dist1.^2))*sqrt(sum(dist2.^2))))*180/π, digits=3), " in degrees" )
    return acos(sum(dist1.*dist2)/(sqrt(sum(dist1.^2))*sqrt(sum(dist2.^2))))
end

"""
$(TYPEDSIGNATURES)

Returns the distance between two atoms given in cartesian coordinates. 

"""
function getdist(centeratom::Vector{<:Real}, branch::Vector{<:Real})
    #Take dot product of the two displacement vectors and divide by the lengths 
    #Then take arccos of answer
    dist = sqrt(sum((branch-centeratom).^2))
    return dist
end

"""
$(TYPEDSIGNATURES)

Returns the distance between two ions read from JDFTX formatted files

`ionpos` : Name of JDFTX ionic positions file

`lattice` : Name of JDFTX lattice file

`numc` : Index of first ion

`numb` : Index of second ion 

Note that the distances are likely returned in Bohr since the the JDFTX lattice is likely written in atomic units
"""
function getdist(ionpos::String, lattice::String, numc::Integer, numb::Integer)
    cartesiancoords = tocartesians(ionpos, lattice)
    return getdist(cartesiancoords[numc], cartesiancoords[numb])
end

"""
$(TYPEDSIGNATURES)

Returns the angle between `v1` and `v2` where `v1` is the vector connecting `numc` and `numb1`
and `v2` is the vector connecting `numc` and `numb2`

`ionpos` : Name of JDFTX formatted ionic positions file

`lattice` : Name of JDFTX formatted lattice file 

`numc` : Index of the center atom

`numb1` : Index of the first branch atom

`numb2` : Index of the second branch atom

"""
function getangle(ionpos::String, lattice::String, numc::Integer, numb1::Integer, numb2::Integer)
    cartesiancoords = tocartesians(ionpos, lattice)
    getangle(cartesiancoords[numc], cartesiancoords[numb1], cartesiancoords[numb2])
end

"""
$(TYPEDSIGNATURES)
Returns the ionic position provided in lattice coordinates in cartesian coordinates. 

`lattype` denotes the lattice type (lattice vectors being columns or rows of provided 3x3 matrix)- with jdftx being the default. 

"""
function tocartesian(ion::Vector{<:Real}, lattice::Array{<:Real, 2}; lattype::Symbol=:jdftx)
    if lattype == :jdftx
        println("lattype interpreted as jdftx convention (columns are lattice vectors)")
        return lattice*ion
    elseif lattype == :normal
        println("Lattype interpreted as normal (rows are lattice vectors) ")
        return transpose(lattice)*ion
    else 
        error("lattype is unrecognized: must be one of :jdftx or :normal")
    end
end

"""
$(TYPEDSIGNATURES)
Returns the ionic positions provided in lattice coordinates in cartesian coordinates. 

`lattype` denotes the lattice type (lattice vectors being columns or rows of provided 3x3 matrix)- with jdftx being the default. 
"""
function tocartesians(ions::Vector{<:Vector{<:Real}}, lattice::Array{<:Real, 2}; lattype::Symbol=:jdftx)
    cartions = Vector{Vector{Float64}}()
    for ion in ions
        push!(cartions, tocartesian(ion, lattice; lattype))
    end
    return cartions
end

"""
$(TYPEDSIGNATURES)
Returns the provided ionic positions (given in the form of an ionic positions file by JDFTX conventions)
in cartesian coordinates. The lattice file is also read from a lattice file assumed to be written in JDFTX conventions. 

`ionpos` : A JDFTX ionic positions file

`lattice` : A JDFTX lattice file
"""
function tocartesians(ionpos::String, lattice::String)
    readions = Vector{Vector{Float64}}()
    latt = Array{Float64, 2}(undef, (3, 3))
    for line in readlines(ionpos)
        push!(readions, parse.(Float64, String.(split(line)[3:5])))
    end
    counter = 1
    for line in readlines(lattice)
        try
            latt[counter, :] = parse.(Float64, String.(split(replace(line, "\\"=>""))[1:3])) #Parses the lattice file
            counter += 1
        catch
        end
    end
    println("Lattice is: ", latt)
    return tocartesians(readions, latt; lattype=:jdftx)
end