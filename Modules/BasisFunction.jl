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

    # gaussian basis function
    Ψ = exp.(-shape_parameter^2 .* dist .^ 2)

    # normalization - eps() avoid division by zero

    return Ψ ./ (sum(Ψ; dims=2) .+ eps())
end

# ==============================================================================
end