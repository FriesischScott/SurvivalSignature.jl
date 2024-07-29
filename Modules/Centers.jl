module Centers

# ==============================================================================

using IterTools

# ==============================================================================

using ..Structures: GridCenters

# ==============================================================================

export generateCenters

# ============================== METHODS =======================================

function generateCenters(method::GridCenters, state_vectors::Matrix, threshold::Float64)
    # grid centers
    lb = minimum(state_vectors; dims=2)
    ub = maximum(state_vectors; dims=2)

    ranges = [range(l, u, c) for (l, u, c) in zip(lb, ub, method.centers_interval)]

    return hcat(
        [
            [c...] for
            c in IterTools.Iterators.product(ranges...) if sum(c .- 1) > threshold
        ]...,
    )
end

# ==============================================================================
end