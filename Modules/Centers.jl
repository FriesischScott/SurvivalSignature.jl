module Centers

# ==============================================================================

using IterTools

# ==============================================================================

using ..Structures: GridCenters

# ==============================================================================

export generateCenters

# ============================== METHODS =======================================

function generateCenters(
    method::String, state_vectors::AbstractArray, threshold::Number, ci::AbstractArray
)
    method = lowercase(method)
    mode = if method == "grid-aligned"
        GridCenters()
    else
        error("Unrecognized Method: $method")
    end

    return generateCenters(mode, state_vectors, threshold, ci)
end

# ==============================================================================

function generateCenters(
    method::GridCenters, state_vectors::AbstractArray, threshold::Number, ci::AbstractArray
)
    # grid centers
    lb = minimum(state_vectors; dims=2)
    ub = maximum(state_vectors; dims=2)

    ranges = [range(l, u, c) for (l, u, c) in zip(lb, ub, ci)]

    return hcat(
        [
            [c...] for
            c in IterTools.Iterators.product(ranges...) if sum(c .- 1) > threshold
        ]...,
    ),
    method
end

# ==============================================================================
end