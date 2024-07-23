module StartingPoints

# ==============================================================================

using NearestNeighbors
using IterTools

# ==============================================================================
#include("Structures.jl")
using ..Structures: Points
using ..Structures: GridStart

#include("SurvivalSignatureUtils.jl")
using ..SurvivalSignatureUtils
# ==============================================================================

export generateStartingPoints

# ============================== METHODS =======================================

function generateStartingPoints(method::String, state_vectors, types)
    method = lowercase(method)
    mode = if method == "grid-aligned"
        GridStart()
    else
        error("Unrecognized Method: $method")
    end

    return generateStartingPoints(mode, state_vectors, types)
end

function generateStartingPoints(method::GridStart, state_vectors, types)
    # grid startng points
    lb = minimum(state_vectors; dims=2)
    ub = maximum(state_vectors; dims=2)

    tree = NearestNeighbors.KDTree(state_vectors)

    ranges = [range(0.0, 1.0; length=l) for l in fill(5, length(types))]  # creates range objects

    starting_points = mapreduce(t -> [t...], hcat, IterTools.Iterators.product(ranges...))

    starting_points = (starting_points .* (ub .- lb) .+ lb)   # scales Xn based on ub and lb

    idx, _ = NearestNeighbors.nn(tree, starting_points)  # nearest neighbor (nn)
    idx = unique(idx)

    # in case two have the same nearest neighbor
    # finds nearest neighbors to the grid aligned starting points - use the C values instead.

    return Points(state_vectors[:, idx], idx, nothing, nothing), method
end
# ==============================================================================

end