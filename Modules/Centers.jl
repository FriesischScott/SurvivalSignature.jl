module Centers

# ==============================================================================

using IterTools

# ==============================================================================
export generateCenters
# ============================== METHODS ======================================

function gridAlignedCenters(
    state_vectors::AbstractArray, threshold::Number, ci::AbstractArray
)
    lb = minimum(state_vectors; dims=2)
    ub = maximum(state_vectors; dims=2)

    ranges = [range(l, u, c) for (l, u, c) in zip(lb, ub, ci)]

    return hcat(
        [
            [c...] for
            c in IterTools.Iterators.product(ranges...) if sum(c .- 1) > threshold
        ]...,
    ),
    ranges
end

# ============================ SELECTION =======================================

function selectCentersMethod(method::String)
    methods_dict = Dict("grid-aligned" => gridAlignedCenters)

    if haskey(methods_dict, lowercase(method))
        methods_dict[lowercase(method)]
    else
        error("Unsupported Method: $method")
    end
end

function generateCenters(
    method::String, state_vectors::AbstractArray, threshold::Number, ci::AbstractVector
)
    func::Function = selectCentersMethod(method)

    return func(state_vectors, threshold, ci)
end

# ==============================================================================
end