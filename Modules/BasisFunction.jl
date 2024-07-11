module BasisFunction

# ==============================================================================

using IterTools
using LinearAlgebra

# ==============================================================================

#include("Structures.jl")
using ..Structures: Points

#include("SurvivalSignatureUtils.jl")
using ..SurvivalSignatureUtils

# ==============================================================================
export basis
# ==============================================================================

# Define a constant dictionary for function types
#   FUNCTION_METHODS = Dict(
#       "gaussian" => gaussian, "matern" => matern, "behrensdorf" => behrensdorf
#   )

function basis(
    function_method::String,
    shape_parameter::Union{Number,AbstractArray},   # behrensdorf (array)
    X::AbstractArray,                               # coordinates
    Y::AbstractArray,                               # centers
    smoothness_factor::Int=1,
)
    Ψ = if function_method == "gaussian"
        # Compute pairwise Euclidean distances
        dist = [
            LinearAlgebra.norm(x .- c) for
            (x, c) in IterTools.product(eachcol(X), eachcol(Y))
        ]

        gaussian(shape_parameter, dist)
    elseif function_method == "matern"

        # Compute pairwise Euclidean distances
        dist = [
            LinearAlgebra.norm(x .- c) for
            (x, c) in IterTools.product(eachcol(X), eachcol(Y))
        ]

        matern(shape_parameter, dist, smoothness_factor)
    elseif function_method == "behrensdorf"
        behrensdorf(shape_parameter, X, Y)
    else
        error("Unsupported function type: $function_type")
    end

    return Ψ ./ (sum(Ψ; dims=2) .+ eps()) # normalization - eps() avoid division by zero
end

function gaussian(shape_parameter::Union{Number,AbstractArray}, X::AbstractMatrix)
    return exp.(-shape_parameter^2 .* X .^ 2)
end

function matern(
    shape_parameter::Union{Number,AbstractArray}, X::AbstractMatrix, smoothness_factor::Int
)
    if smoothness_factor == 1
        return exp.(-shape_parameter .* X)
    elseif smoothness_factor == 2
        return (1 .+ shape_parameter .* X) .* exp.(-shape_parameter .* X)
    elseif smoothness_factor == 3
        return (
            (3 .+ 3 .* shape_parameter .* X + shape_parameter^2 .* X .^ 2) .*
            exp.(-shape_parameter .* X)
        )
    else
        error("Smoothness Factor MUST BE 1, 2, or 3")
    end
end

function behrensdorf(
    shape_parameter::Union{Number,AbstractArray}, X::AbstractArray, centers::AbstractArray
)
    Ψ =
        exp.(
            -[
                sum((x .- c) .^ 2 ./ 2 .* shape_parameter .^ 2) for
                (x, c) in Iterators.product(eachcol(X), eachcol(centers))
            ]
        )
    return Ψ
end

# test behrensdorf and gaussian
# ==============================================================================
end