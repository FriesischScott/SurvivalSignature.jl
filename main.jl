using SurvivalSignature
using ProgressMeter
using JLD2              # needed for load
# ==============================================================================

include("Modules/Import.jl")
using .Import

using ..SurvivalSignatureUtils

using ..Structures
using ..Structures: System, Simulation, Model, Methods
using ..Structures: SystemMethod, GridSystem
using ..Structures: SimulationType, MonteCarloSimulation, IntervalPredictorSimulation
using ..Structures: BasisFunctionMethod, Gaussian
using ..Structures: ShapeParameterMethod, Hardy, Rippa, DirectAMLS, IndirectAMLS
using ..Structures: CentersMethod, GridCenters
using ..Structures: StartingMethod, GridStart
using ..Structures: ErrorType, RMSE, RAE, NORM, NRMSE
using ..Structures: AdaptiveRefinementMethod, TEAD

using ..StructureCompilation
using ..Error: calculateError
using ..Systems
using ..Visualization: plotError

function main()
    # ================================ INPUTS ==================================
    percolation::Bool = true
    verbose::Bool = true        # used to turn on and off print statements during 'simulate'

    # [ GridSystem() ]
    system_type::SystemMethod = GridSystem((15, 15))

    # Simulation Parameters
    samples::Union{Int,Vector{Int}} = [1, 10, 25, 50, 100, 250, 500, 1000]
    covtol::Float64 = 1e-3                       # coeficient of varriation tolerance
    wtol::Float64 = 1e-4                         # weight change tolerance

    ci::Vector{Int} = [15, 15]               # centers interval - dims must match 
    #                                            # number of types

    # METHODS
    # -------------------------------------------------------------------------
    # -------------------------------------------------------------------------
    # [MonteCarloSimulation(), IntervalPredictorSimulation()]
    simulation_method::SimulationType = IntervalPredictorSimulation()
    # [ GridStart() ]
    starting_points_method::StartingMethod = GridStart()
    # [ GridCenters() ]
    centers_method::CentersMethod = GridCenters(ci)
    # [ Norm() ]
    weight_change_method::ErrorType = NORM()
    # [ Hardy(), Rippa(), DirectAMLS() ] # can be a Vector
    shape_parameter_method::Union{Vector{ShapeParameterMethod},ShapeParameterMethod} = [
        Hardy(), Rippa(), DirectAMLS(), IndirectAMLS()
    ]

    # [ Gaussian() ] 
    basis_function_method::Union{Vector{BasisFunctionMethod},BasisFunctionMethod} = Gaussian()
    # [ TEAD() ]
    adaptive_refinement_method::AdaptiveRefinementMethod = TEAD()

    # Error
    error_type::Union{Vector{ErrorType},ErrorType} = [RMSE(), RAE()]

    # ======================== STRUCT REFINEMENT ===============================

    sys::System = Systems.generateSystem(system_type; percolation_bool=percolation)

    sims::Union{Vector{Simulation},Simulation} = StructureCompilation.compileSimulation(
        samples, covtol, wtol
    )

    methods::Union{Vector{Methods},Methods} = StructureCompilation.compileMethods(
        simulation_method,
        starting_points_method,
        centers_method,
        weight_change_method,
        shape_parameter_method,
        basis_function_method,
        adaptive_refinement_method,
    )

    # =============================== SIMULATE =================================

    signatures::Union{Matrix{Model},Vector{Model},Model} = Simulate.simulate(
        methods, sys, sims; verbose=verbose
    )

    # ============================= "TRUE SOLUTION" ============================
    # sim_mc = Structures.Simulation(1000, covtol, wtol)
    # method_mc = Structures.Methods(
    #     MonteCarloSimulation(), nothing, nothing, nothing, nothing, nothing, nothing
    # )
    # signature_mc = Simulate.simulate(method_mc, sys, sim_mc; verbose=verbose)

    #true values - apparently this outputs a Φ - which is the true value solutions
    @load "demo/data/grid-network-15x15-MC-10000.jld2"

    # =============================== ERROR ====================================
    println("Calculating Error...")

    # signatures, errors = Error.calculateError(
    #     error_type, signatures, signature_mc.Phi.solution; verbose=false
    # )

    signatures, errors = Error.calculateError(error_type, signatures, Φ; verbose=false)

    SurvivalSignatureUtils._print(errors)

    println("Errors Calculated.\n")

    plt = Visualization.plotError(signatures, errors)

    display(plt)

    # Timings:

    println("Time (s):")

    metrics = Matrix{Metrics}(undef, size(signatures))
    for (i, signature) in enumerate(signatures)
        metrics[i] = signature.metrics
    end

    values = Matrix{Float64}(undef, size(metrics))
    for (i, metric) in enumerate(metrics)
        values[i] = metric.time
    end

    SurvivalSignatureUtils._print(values)

    # SurvivalSignatureUtils._print(signatures[1, end].Phi.solution)
    # SurvivalSignatureUtils._print(signatures[2, end].Phi.solution)
    # SurvivalSignatureUtils._print(signatures[3, end].Phi.solution)

    return nothing
end
# =============================== RUN ==========================================
main()
