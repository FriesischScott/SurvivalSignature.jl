module Structures

# ==============================================================================
using ..SurvivalSignatureUtils

# ==============================================================================
export System, Simulation, StringMethods, Methods
export Model, PredictorModel
export Gaussian, Matern, MultiQuadaratic, BehrensdorfBasis
export MonteCarloSimulation, RadialBasisSimulation, IntervalPredictorSimulation
export Hardy, Franke, Kuo, Rippa, BehrensdorfShape
export GridCenters
export GridStart

# ==============================================================================
struct System
    adj::AbstractArray
    connectivity::Any
    types::Dict{Int64,Vector{Int64}}
    percolation::Bool

    # Constructor with default values 
    function System(
        adj::AbstractArray,
        connectivity::Any,
        types::Dict{Int64,Vector{Int64}},
        percolation::Bool=true,
    )
        return new(adj, connectivity, types, percolation)
    end
end

mutable struct Simulation
    samples::Int
    variation_tolerance::Float64
    weight_change_tolerance::Float64
    confidence_interval::Vector{Int}
    threshold::Union{Number,Nothing}  # nothing prior to percolation

    # Constructor with default values
    function Simulation(
        samples::Int=1000,
        variation_tolerance::Float64=1e-3,
        weight_change_tolerance::Float64=1e-3,
        confidence_interval::Vector{Int}=[15, 15],
        threshold::Union{Number,Nothing}=nothing,
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

@allow_nothing struct StringMethods
    simulation_method::String
    starting_points_method::String
    centers_method::String
    weight_change_method::String
    shape_parameter_method::String
    basis_function_method::String
    smoothness_factor::Int

    # Constructor with default values
    function StringMethods(
        simulation_method::String="monte-carlo",
        starting_points_method::String="grid",
        centers_method::String="grid",
        weight_change_method::String="norm",
        shape_parameter_method::String="hardy",
        basis_function_method::String="gaussian",
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

# Adding comments within a struct with the macro @allow_nothing causes issues.
@allow_nothing struct Points
    coordinates::Array
    idx::Union{Number,Vector,Matrix}
    solution::Union{Number,Vector,Matrix}
    confidence::Union{Number,Vector,Matrix}
end

# Coordinates will almost always be an Array, but Any allows for variable types 
# used in some modules where it is unreasonable to predict or write them all 

# ============================== BASIS FUNCTIONS ===============================

struct Gaussian end

struct Matern
    smoothness_factor::Int
end

struct MultiQuadaratic
    degree::Int
end

struct BehrensdorfBasis end

# ============================== SIMULATIONS ===================================

struct MonteCarloSimulation end

struct RadialBasisSimulation end

struct IntervalPredictorSimulation end

# ============================= SHAPE PARAMETERS ================================

struct Hardy
    points::Array
end

struct Franke
    points::Array
end

struct Kuo
    points::Array
end

struct Rippa
    starting_points::Points
    centers::Array
end

struct BehrensdorfShape
    confidence_interval::Vector{Int}
    upper::Array
    lower::Array
end

# ================================ CENTERS =====================================

struct GridCenters end

# ================================ STARTING ====================================

struct GridStart end

# ==============================================================================

struct Methods
    simulation_method::Union{
        MonteCarloSimulation,RadialBasisSimulation,IntervalPredictorSimulation
    }
    starting_points_method::Union{GridStart,Nothing}
    centers_method::Union{GridCenters,Nothing}
    weight_change_method::Union{String,Nothing}
    shape_parameter_method::Union{Hardy,Franke,Kuo,Rippa,BehrensdorfShape,Nothing}
    basis_function_method::Union{Gaussian,Matern,MultiQuadaratic,BehrensdorfBasis,Nothing}
end

# =========================== MODELS ===========================================
# Used for Interval Predictor Models
struct PredictorModel
    evaluated_points::Points
    centers::Array
    shape_parameter::Union{Number,Array}
    weights::Array
    w_u::Vector{Float64}
    w_l::Vector{Float64}
end

mutable struct Model
    Phi::Points
    model::Union{PredictorModel,Nothing}                # Specific to the type of simulation
    sys::System                                 # For post-processing access
    sim::Simulation
    method::Methods
    str_method::StringMethods
    errors::Union{Nothing,Dict{String,Float64}}
end

# ==============================================================================

end
