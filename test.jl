
using Random

include("Modules/Import.jl")
using .Import

using ..SurvivalSignatureUtils: _print

function rand_dict()
    rand1 = Random.rand()
    rand2 = Random.rand()

    return Dict("error1" => rand1, "error2" => rand2, "error3" => "works")
end

# Initialize the matrix with the correct type
matrix1 = Matrix{Dict{String,Union{String,Float64}}}(undef, 25, 25)
#matrix1 = Matrix{Int}(undef, 5, 5)

# Assign a dictionary to each index of the matrix
for i in 1:size(matrix1, 1)
    for j in 1:size(matrix1, 2)
        matrix1[i, j] = rand_dict()
    end
end

SurvivalSignatureUtils._print(matrix1)