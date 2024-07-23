module BasisFunction

# ==============================================================================

using IterTools
using LinearAlgebra

# ==============================================================================

#include("Structures.jl")
using ..Structures: Points
using ..Structures: StringMethods
using ..Structures: Gaussian, Matern, MultiQuadaratic, BehrensdorfBasis

#include("SurvivalSignatureUtils.jl")
using ..SurvivalSignatureUtils

# ==============================================================================
export basis, selectBasisFunctionMethod
# ==============================================================================

function basis(
    method::StringMethods,
    shape_parameter::Union{Number,AbstractVector},
    coordinates::AbstractArray,
    centers::AbstractArray,
)
    method = lowercase(method.basis_function_method)
    mode = if method == "gaussian"
        Gaussian()
    elseif method == "behrensdorf"
        BehrensdorfBasis()
    elseif method == "matern"
        Matern(method.smoothness_factor)
    elseif method == "mq"
        MultiQuadaratic(method.smoothness_factor)
    end

    return basis(mode, shape_parameter, coordinates, centers)
end

# ==============================================================================

function basis(
    method::Gaussian,
    shape_parameter::Number,
    coordinates::AbstractArray,
    centers::AbstractArray,
)
    dist = [
        LinearAlgebra.norm(x .- c) for
        (x, c) in IterTools.product(eachcol(coordinates), eachcol(centers))
    ]

    # gaussian basis function
    Ψ = exp.(-shape_parameter^2 .* dist .^ 2)

    # normalization - eps() avoid division by zero

    return Ψ ./ (sum(Ψ; dims=2) .+ eps()), method
end
function basis(
    method::MultiQuadaratic,
    shape_parameter::Number,
    coordinates::AbstractArray,
    centers::AbstractArray,
)
    # mutliquadratic basis functions
    dist = [
        LinearAlgebra.norm(x .- c) for
        (x, c) in IterTools.product(eachcol(coordinates), eachcol(centers))
    ]

    Ψ = if method.degree == 1
        sqrt.(dist .^ 2 .+ shape_parameter^2)
    elseif method.degree == 2
        sqrt.(dist .^ 2 .* shape_parameter^2 .+ 1)
    else
        error("Degree must be 1 or 2 - Inputed: $(method.degree)")
    end
    # normalization - eps() avoid division by zero
    return Ψ ./ (sum(Ψ; dims=2) .+ eps()), method
end

function basis(
    method::Matern,
    shape_parameter::Number,
    coordinates::AbstractArray,
    centers::AbstractArray,
)
    # matern basis function
    dist = [
        LinearAlgebra.norm(x .- c) for
        (x, c) in IterTools.product(eachcol(coordinates), eachcol(centers))
    ]

    Ψ = if method.smoothness_factor == 1
        exp.(-shape_parameter .* dist)
    elseif method.smoothness_factor == 2
        (1 .+ shape_parameter .* dist) .* exp.(-shape_parameter .* dist)
    elseif method.smoothness_factor == 3
        (3 .+ 3 .* shape_parameter .* dist + shape_parameter^2 .* dist .^ 2) .*
        exp.(-shape_parameter .* dist)
    else
        error("Smoothness Factor MUST BE 1, 2, or 3 - Inputed: $(method.smoothness_factor)")
    end
    # normalization - eps() avoid division by zero
    return Ψ ./ (sum(Ψ; dims=2) .+ eps()), method
end

function basis(
    method::BehrensdorfBasis,
    shape_parameter::Union{Number,AbstractVector{Float64}},
    coordinates::AbstractArray,
    centers::AbstractArray,
)
    # behrensdorf basis function - should be the same as gaussian - testing necessary
    Ψ =
        exp.(
            -[
                sum((x .- c) .^ 2 ./ 2shape_parameter .^ 2) for
                (x, c) in Iterators.product(eachcol(coordinates), eachcol(centers))
            ]
        )

    # normalization - eps() avoid division by zero
    return Ψ ./ (sum(Ψ; dims=2) .+ eps()), method
end

# ==============================================================================
end