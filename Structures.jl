module Structures

# ==============================================================================

#include("SurvivalSignatureUtils.jl")
using ..SurvivalSignatureUtils

# ==============================================================================

export System, Simulation, Method, Points, PredictorModel

# ==============================================================================

struct System
    adj::AbstractArray
    connectivity::Any
    types::Dict
    percolation::Bool

    # Constructor with default values 
    function System(
        adj::AbstractArray, connectivity::Any, types::Dict, percolation::Bool=true
    )
        return new(adj, connectivity, types, percolation)
    end
end

mutable struct Simulation # mutable to add threshold
    samples::Int
    variation_tolerance::Float64
    weight_change_tolerance::Float64
    confidence_interval::AbstractVector{Int}
    threshold::Union{Nothing,Number}  # nothing prior to percolation

    # Constructor with default values
    function Simulation(
        samples::Int=1000,
        variation_tolerance::Float64=1e-3,
        weight_change_tolerance::Float64=1e-3,
        confidence_interval::AbstractVector{Int}=[15, 15],
        threshold::Union{Nothing,Number}=nothing,
    )
        return new(
            samples,
            variation_tolerance,
            weight_change_tolerance,
            confidence_interval,
            threshold,
        )
    end
end

mutable struct Method
    simulation_method::String
    starting_points_method::String
    centers_method::String
    weight_change_method::String
    shape_parameter_method::String
    basis_function_method::String
    smoothness_factor::Int

    # Constructor with default values
    function Method(
        simulation_method::String="monte-carlo",
        starting_points_method::String="grid-aligned",
        centers_method::String="grid-aligned",
        weight_change_method::String="norm",
        shape_parameter_method::String="behrensdorf",
        basis_function_method::String="behrensdorf",
        smoothness_factor::Int=1,
    )
        return new(
            simulation_method,
            starting_points_method,
            centers_method,
            weight_change_method,
            shape_parameter_method,
            basis_function_method,
            smoothness_factor,
        )
    end
end

# adding comments within a struct with the macro @allow_nothing causes issues.
@allow_nothing struct Points
    coordinates::Any
    idx::Union{Number,AbstractArray}
    solution::Union{Number,AbstractArray}
    confidence::Union{Number,AbstractArray}
end

# coordinates will almost always be an Array, but any allows for varible types 
# used in some modules where it is unreasonable to predict or write them all 

# =========================== MODELS ===========================================
struct Model
    Phi::Points
    model::Any                # specific to the type of simulation
    sys::System               # for post-proccessing access
    sim::Simulation
    method::Method
end

# used for Interval Predictor Models
struct PredictorModel
    evaluated_points::Points
    centers::AbstractArray
    shape_parameter::Union{Number,AbstractArray}
    weights::AbstractArray
    w_u::Vector             # not sure what 'u' stands for
    w_l::Vector             # not sure what 'l' stands for
end
# ==============================================================================

end