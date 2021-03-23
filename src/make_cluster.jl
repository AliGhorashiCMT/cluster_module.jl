
function make_xsf()
run(`createXSF`)
end

function make_test_xsf()
path_to_test_si = joinpath(@__DIR__, "../examples/si.out")
run(`cat $path_to_test_si`)
run(`cp $path_to_test_si si-copy.out`)
run(`createXSF si-copy.out si.xsf `)
run(`open si.xsf`)
end