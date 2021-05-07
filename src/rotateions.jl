"""
$(TYPEDSIGNATURES)
Rotate an ion with respect to a given coordinate. Supply this function with a position 
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
Returns the perpendicular vector to the plane created by three atoms
"""
function findperp(v1::Vector{<:Real}, v2::Vector{<:Real}, v3::Vector{<:Real})
    w1 = (v2-v1)
    w2 = (v3-v1)
    cross(w1, w2)
end

"""
$(TYPEDSIGNATURES)
"""
function extendbond(extendfrom::Vector{<:Real}, extendto::Vector{<:Real}, newlength::Real, ions::Vector{<:Vector{<:Real}})
    vecextend = extendto - extendfrom
    vecextend *= newlength./sqrt(sum(vecextend.^2))
    return [vecextend, [vecextend+extendfrom-extendto+ion for ion in ions]...]
end

"""
$(TYPEDSIGNATURES)
"""
function extendbond(extendfrom::Vector{<:Real}, extendto::Vector{<:Real}, newlength::Real)
    vecextend = extendto - extendfrom
    vecextend *= newlength./sqrt(sum(vecextend.^2))
    return vecextend+extendto
end

"""
$(TYPEDSIGNATURES)
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
For rotating about z axis 
"""
function rotateinz(pos::Vector{<:Real}, θ::Real)
    round.([cos(θ) sin(θ) 0; -sin(θ) cos(θ) 0; 0 0 1]*pos, digits=4)
end