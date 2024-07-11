module Systems

# ==============================================================================
using LinearAlgebra
# ==============================================================================
using ..Structures: System

# needed for grid-network
include("../src/util.jl")
# needed for s_t_connectivity
include("../src/structurefunctions.jl")
# ==============================================================================

export generateSystem

# ==============================================================================

function selectSystem(method::String)
    methods_dict = Dict("grid" => gridSystem)

    if haskey(methods_dict, lowercase(method))
        return methods_dict[lowercase(method)]
    else
        error("Unsupported Method: $method")
    end
end

function generateSystem(method::String, dims...; percolation_bool::Bool=true)
    func = selectSystem(method)

    adj, connectivity, types = func(dims...)
    return System(adj, connectivity, types, percolation_bool)
end

# ==============================================================================
function gridSystem(dims...)
    adj = gridnetwork(dims...)
    φ = s_t_connectivity([1:prod(dims);], [1], [prod(dims)])
    types = Dict(1 => collect(1:2:prod(dims)), 2 => collect(2:2:prod(dims))) # alternate every other type

    return adj, φ, types
end

# ==============================================================================
end
