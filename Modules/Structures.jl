module Structures

# ==============================================================================

using ..SurvivalSignatureUtils

# ==============================================================================
export System, Simulation, Methods

export SimulationModel, MonteCarloModel, PredictorModel, Model
export SimulationType, MonteCarloSimulation, IntervalPredictorSimulation

export BasisFunctionMethod, Gaussian
export ShapeParameterMethod, Hardy, Rippa, DirectAMLS, IndirectAMLS
export CentersMethod, GridCenters
export StartingMethod, GridStart
export SystemMethod, GridSystem
export AdaptiveRefinementMethod, TEAD
export Metrics

export ErrorType, RMSE, RAE, NORM, NRMSE

# ==============================================================================
struct System
    adj::Matrix
    connectivity::Any
    types::Dict{Int64,Vector{Int64}}
    percolation::Bool

    # Constructor with default values 
    function System(
        adj::Matrix, connectivity::Any, types::Dict{Int,Vector{Int}}, percolation::Bool=true
    )
        return new(adj, connectivity, types, percolation)
    end
end

mutable struct Simulation
    samples::Int
    variation_tolerance::Float64
    weight_change_tolerance::Float64
    threshold::Union{Float64,Nothing}  # nothing prior to percolation

    # Constructor with default values
    function Simulation(
        samples::Int=1000,
        variation_tolerance::Float64=1e-3,
        weight_change_tolerance::Float64=1e-3,
        threshold::Union{Float64,Nothing}=nothing,
    )
        return new(samples, variation_tolerance, weight_change_tolerance, threshold)
    end
end

mutable struct Points
    coordinates::Union{Array,Nothing}
    idx::Union{Number,Vector,Matrix,Nothing}
    solution::Union{Float64,Vector,Matrix,Nothing}
    confidence::Union{Float64,Vector,Matrix,Nothing}
end

# ============================== BASIS FUNCTIONS ===============================

abstract type BasisFunctionMethod end

struct Gaussian <: BasisFunctionMethod end

# ============================== SIMULATIONS ===================================

abstract type SimulationType end

struct MonteCarloSimulation <: SimulationType end

struct IntervalPredictorSimulation <: SimulationType end

# ============================= SHAPE PARAMETERS ================================

abstract type ShapeParameterMethod end
struct Hardy <: ShapeParameterMethod end
struct Rippa <: ShapeParameterMethod end
struct DirectAMLS <: ShapeParameterMethod
    max_iterations::Int
    tolerance::Float64

    function DirectAMLS(
        max_iterations::Union{Nothing,Int}=10, tolerance::Union{Nothing,Float64}=1e-6
    )
        return new(max_iterations, tolerance)
    end
end
struct IndirectAMLS <: ShapeParameterMethod
    max_iterations::Int
    tolerance::Float64

    function IndirectAMLS(
        max_iterations::Union{Nothing,Int}=10, tolerance::Union{Nothing,Float64}=1e-6
    )
        return new(max_iterations, tolerance)
    end
end

# struct Franke <: ShapeParameterMethod end   # only proven with MQ 
# struct Kuo <: ShapeParameterMethod end      # only proven with MQ

# ================================ SYSTEMS =====================================

abstract type SystemMethod end

struct GridSystem <: SystemMethod
    dims::Tuple
end

# ================================ SYSTEMS =====================================

abstract type AdaptiveRefinementMethod end

struct TEAD <: AdaptiveRefinementMethod end

# ================================ CENTERS =====================================

abstract type CentersMethod end

struct GridCenters <: CentersMethod
    centers_interval::Vector{Int}       # number of centers in each dimension

    function GridCenters(centers_interval::Vector{Int}=[15, 15])
        return new(centers_interval)
    end
end

# ================================ CENTERS =====================================

abstract type ErrorType end

struct RMSE <: ErrorType end
struct RAE <: ErrorType end
struct NORM <: ErrorType end
struct NRMSE <: ErrorType end

# ================================ STARTING ====================================

abstract type StartingMethod end

struct GridStart <: StartingMethod end

# ==============================================================================

struct Methods
    simulation_method::Union{SimulationType,Nothing}
    starting_points_method::Union{StartingMethod,Nothing}
    centers_method::Union{CentersMethod,Nothing}
    weight_change_method::Union{ErrorType,Nothing}
    shape_parameter_method::Union{ShapeParameterMethod,Nothing}
    basis_function_method::Union{BasisFunctionMethod,Nothing}
    adaptive_refinement_method::Union{AdaptiveRefinementMethod,Nothing}
end
# ================================ DATA ========================================
struct Metrics
    time::Float64
    # model_evaluations
end

# =========================== MODELS ===========================================
abstract type SimulationModel end

struct MonteCarloModel <: SimulationModel end

struct PredictorModel <: SimulationModel
    evaluated_points::Points
    centers::Array
    shape_parameter::Float64
    weights::Array
    w_u::Vector{Float64}
    w_l::Vector{Float64}
end

mutable struct Model
    Phi::Points
    model::SimulationModel                          # Specific to the type of simulation
    sys::System                                     # For post-processing access
    sim::Simulation
    method::Methods
    metrics::Union{Metrics,Nothing}
    errors::Union{Dict{String,Float64},Nothing}     # starts nothing, then populated
end

# ==============================================================================

end
