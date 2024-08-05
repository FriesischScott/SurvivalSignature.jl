module BasisFunction

# ==============================================================================

using IterTools
using LinearAlgebra

# ==============================================================================

using ..Structures: Gaussian
using ..SurvivalSignatureUtils

# ==============================================================================
export basis
# ==============================================================================

function basis(
    method::Gaussian,
    shape_parameter::Float64,
    coordinates::Union{Matrix,Vector},
    centers::Matrix,
)
    dist = [
        LinearAlgebra.norm(x .- c) for
        (x, c) in IterTools.product(eachcol(coordinates), eachcol(centers))
    ]

    Ψ = exp.(-dist .^ 2 ./ (2 * shape_parameter^2))

    # normalization - eps() avoid division by zero
    return Ψ ./ (sum(Ψ; dims=2) .+ eps())
end

function basis(method::Gaussian, shape_parameter::Float64, dist::Float64)
    # used in indirect AMLS shape parameter method

    Ψ = exp.(-dist .^ 2 ./ (2 * shape_parameter^2))

    # normalization not necessary since Ψ is scalar
    return Ψ
end

end

# ==============================================================================
