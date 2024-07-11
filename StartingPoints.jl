module StartingPoints

# ==============================================================================

using NearestNeighbors
using IterTools

# ==============================================================================
#include("Structures.jl")
using ..Structures: Points

#include("SurvivalSignatureUtils.jl")
using ..SurvivalSignatureUtils

# ============================== METHODS =======================================
function gridAlignedStartingPoints(state_vectors, types)
    lb = minimum(state_vectors; dims=2)
    ub = maximum(state_vectors; dims=2)

    tree = NearestNeighbors.KDTree(state_vectors)

    ranges = [range(0.0, 1.0; length=l) for l in fill(5, length(types))]  # creates range objects

    starting_points = mapreduce(t -> [t...], hcat, IterTools.Iterators.product(ranges...))

    starting_points = (starting_points .* (ub .- lb) .+ lb)   # scales Xn based on ub and lb

    idx, _ = nn(tree, starting_points)  # nearest neighbor (nn)
    idx = unique(idx)

    # in case two have the same nearest neighbor
    # finds nearest neighbors to the grid aligned starting points - use the C values instead.

    return Points(state_vectors[:, idx], idx, nothing, nothing)
end

# ============================ SELECTION =======================================

function selectStartingPointsMethod(method::String)
    methods_dict = Dict("grid-aligned" => gridAlignedStartingPoints)

    if haskey(methods_dict, lowercase(method))
        methods_dict[lowercase(method)]
    else
        error("Unsupported Method: $method")
    end
end

function generateStartingPoints(method::String, Ω::AbstractArray, types::Dict)
    func::Function = selectStartingPointsMethod(method)

    return func(Ω, types)
end

# ==============================================================================

export *

end