"""
$(TYPEDSIGNATURES)
"""
function getangle(centeratom::Vector{<:Real}, branch1::Vector{<:Real}, branch2::Vector{<:Real})
    #Take dot product of the two displacement vectors and divide by the lengths 
    #Then take arccos of answer
    dist1 = branch1-centeratom
    dist2 = branch2-centeratom
    println(round(acos(sum(dist1.*dist2)/(sqrt(sum(dist1.^2))*sqrt(sum(dist2.^2))))/π, digits=3), " in units of pi" )
    return acos(sum(dist1.*dist2)/(sqrt(sum(dist1.^2))*sqrt(sum(dist2.^2))))
end

