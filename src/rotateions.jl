"""
$(TYPEDSIGNATURES)
Rotate an ion with respect to a given coordinate. Supply this function with a position 

`pos` : The vector that is to be rotated 

`perp` : The vector about which to rotate. Basically the direction that is unaffected by the rotation.

`θ` : The angle by which to rotate. Note that the angle should be given in radians. 
"""
function rotateion(pos::Vector{<:Real}, perp::Vector{<:Real}, θ::Real)
    length(perp) != 3 && error("Perp must be 3 vector")
    perp ≈ zeros(3) && error("Perp must be nonzero")
    nx, ny, nz = perp/sqrt(sum(perp.^2))
    (nx, ny) == (0, 0) && return(rotateinz(pos, θ))
    eigenvectors = [[nx, ny, nz], [-ny, nx, 0]/sqrt(nx^2+ny^2), [nx*nz, ny*nz, -ny^2-nx^2]/sqrt(nx^2+ny^2)]
    eigenarray = Array{Float64, 2}(undef, (3, 3))
    for (index, col) in enumerate(eachcol(eigenarray))
        eigenarray[:, index]= eigenvectors[index]
    end
    inveigenarray = inv(eigenarray)
    rotatedpos = eigenarray*[1 0 0; 0 cos(θ) sin(θ); 0 -sin(θ) cos(θ)]*inveigenarray*pos
    @assert sum(rotatedpos.^2) ≈ sum(pos.^2) #Check that normalization is left unchanged
    return(round.(rotatedpos, digits=4 ))
    @assert eigenarray*inv(eigenarray) ≈ [1 0 0; 0 1 0; 0 0 1] #Check that inverse is correct
    @assert transpose(eigenarray) ≈ inv(eigenarray) #Check that matrix is orthogonal 
end

"""
$(TYPEDSIGNATURES)
`pos` : Position of ion to be rotated in cartesian coordinates

`aboutwhere` : The effective origin of rotation

`perp` : The direction normal to the rotation

`θ` : The angle by which to rotate

"""
function rotateion(pos::Vector{<:Real}, aboutwhere::Vector{<:Real}, perp::Vector{<:Real}, θ::Real)
    rotateion(pos-aboutwhere, perp, θ) + aboutwhere
end


"""
$(TYPEDSIGNATURES)
Returns the perpendicular vector to the plane created by three atoms. This should give the `perp`
vector to be used by `rotateion`.
"""
function findperp(v1::Vector{<:Real}, v2::Vector{<:Real}, v3::Vector{<:Real})
    w1 = (v2-v1)
    w2 = (v3-v1)
    cross(w1, w2)
end

"""
$(TYPEDSIGNATURES)

Extend (or shorten) a given bond. 

`extendfrom` : The atom position from which the bond should be extended

`extendto` : The atom on the other side of the bond. 

`newlength` : The new bond length

`ions` : The other ions that are connected to the subcluster formed by the bond (Note
        that if a bond is extended all the atoms connected to the bonded atom should be moved as well-
        otherwise their own bonds are changed in direction, length) 
"""
function extendbond(extendfrom::Vector{<:Real}, extendto::Vector{<:Real}, 
    newlength::Real, ions::Vector{<:Vector{<:Real}})

    vecextend = extendto - extendfrom
    vecextend *= newlength./sqrt(sum(vecextend.^2))
    return [vecextend+extendfrom, [vecextend+extendfrom-extendto+ion for ion in ions]...]
end

"""
$(TYPEDSIGNATURES)
"""
function extendbond(extendfrom::Vector{<:Real}, extendto::Vector{<:Real}, newlength::Real)
    vecextend = extendto - extendfrom
    vecextend *= newlength./sqrt(sum(vecextend.^2))
    return vecextend+extendfrom
end

"""
$(TYPEDSIGNATURES)

`poses` : A vector of ionic positions. Note that these positions should be in the cartesian basis- not the
commonly used lattice basis.

`perp` : The direction perpendicular to the rotation. 

`θ` : The angle by which to rotate
"""
function rotateions(poses::Vector{<:Vector{<:Real}}, perp::Vector{<:Real}, θ::Real)
    rotatedposes = Vector{Vector{Float64}}()
    for pos in poses
        push!(rotatedposes, rotateion(pos, perp, θ))
    end
    return rotatedposes
end

"""
$(TYPEDSIGNATURES)

`poses` : A vector of ionic positions. Note that these positions should be in the cartesian basis- not the
commonly used lattice basis.

`aboutwhere` : The origin of the rotation

`perp` : The direction perpendicular to the rotation. 

`θ` : The angle by which to rotate

"""
function rotateions(poses::Vector{<:Vector{<:Real}}, aboutwhere::Vector{<:Real}, perp::Vector{<:Real}, θ::Real)
    rotatedposes = Vector{Vector{Float64}}()
    for pos in poses
        push!(rotatedposes, rotateion(pos, aboutwhere, perp, θ))
    end
    return rotatedposes
end



"""
$(TYPEDSIGNATURES)
For rotating about z axis 

`pos` : The ionic position (in cartesian coordinates) that is to be rotated.

`θ` : The angle by which to rotate about the z axis. 
"""
function rotateinz(pos::Vector{<:Real}, θ::Real)
    round.([cos(θ) sin(θ) 0; -sin(θ) cos(θ) 0; 0 0 1]*pos, digits=4)
end