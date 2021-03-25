
function make_xsf(ionpos::AbstractVector, lattice::Array{<:Real, 2})

run(`createXSF`)

end

function make_test_xsf()
    path_to_test_si = joinpath(@__DIR__, "../examples/si.out");
    print("here")
    if isfile("si-copy.out")
        rm("si-copy.out")
    end
    run(`cp $path_to_test_si si-copy.out`)
    run(`createXSF si-copy.out si.xsf `)
    run(`open si.xsf`)
    rm("si-copy.out")
end