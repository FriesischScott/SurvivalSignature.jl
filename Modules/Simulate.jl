module Simulate

# this module to store the functions relating to simulating a given grid or model.
# ==============================================================================

using ..SurvivalSignatureUtils
using ..Structures: System, Simulation, Methods, Points, Model
using ..Structures: SimulationType, MonteCarloSimulation, IntervalPredictorSimulation
using ..Structures: MonteCarloModel, PredictorModel
using ..Structures: Metrics

using ..Evaluation
using ..StartingPoints
using ..Centers
using ..ShapeParameter
using ..BasisFunction
using ..AdaptiveRefinement
using ..IntervalPredictorModel
using ..Visualization

# convert this to a Module: RadialBasisFunctions (?)
# needed for monotonicity_constraints and lsqr
include("../src/rbf/radialbasisfunctions.jl")

# ==============================================================================
export simulate
# ==============================================================================

function simulate(
    methods::Vector{Methods},
    sys::System,
    sims::Vector{Simulation};
    verbose::Bool=false,
    shape_parameter::Union{Nothing,Float64}=nothing,
)::Matrix{Model}

    # Initialize the signatures array
    signatures = Matrix{Model}(undef, (length(methods), length(sims)))

    for (i, method) in enumerate(methods)   # columns
        for (j, sim) in enumerate(sims) #   # rows
            signatures[i, j] = simulate(
                method.simulation_method,
                sys,
                sim,
                method;
                verbose=verbose,
                shape_parameter=shape_parameter,
            )
        end
    end

    return signatures
end

function simulate(
    methods::Vector{Methods},
    sys::System,
    sim::Simulation;
    verbose::Bool=false,
    shape_parameter::Union{Nothing,Float64}=nothing,
)::Vector{Model}

    # Length check
    len = length(methods)

    #Initialize the signatures array
    signatures = Vector{Model}(undef, len)

    # Iterate over the elements of all vectors simultaneously
    for (i, method) in enumerate(methods)
        signatures[i] = simulate(
            method.simulation_method,
            sys,
            sim,
            method;
            verbose=verbose,
            shape_parameter=shape_parameter,
        )
    end

    return signatures
end

function simulate(
    method::Methods,
    sys::System,
    sims::Vector{Simulation};
    verbose::Bool=false,
    shape_parameter::Union{Nothing,Float64}=nothing,
)::Vector{Model}

    # Initialize the signatures array
    signatures = Vector{Model}(undef, length(sims))

    # Iterate over the elements of all vectors simultaneously
    for (i, sim) in enumerate(sims)
        signatures[i] = simulate(
            method.simulation_method,
            sys,
            sim,
            method;
            verbose=verbose,
            shape_parameter=shape_parameter,
        )
    end

    return signatures
end

function simulate(
    method::Methods,
    sys::System,
    sim::Simulation;
    verbose::Bool=false,
    shape_parameter::Union{Nothing,Float64}=nothing,
)::Model
    return simulate(
        method.simulation_method,
        sys,
        sim,
        method;
        verbose=verbose,
        shape_parameter=shape_parameter,
    )
end

function simulate(
    method::MonteCarloSimulation,
    sys::System,
    sim::Simulation,
    methods::Methods;
    verbose::Bool=false,
    shape_parameter::Union{Nothing,Float64}=nothing,
)::Model
    start_time = time_ns()

    if verbose
        printDetails(sys, sim, methods)
    end
    state_vectors, percolated_state_vectors, sim.threshold = Evaluation.generateStateVectors(
        sys
    )

    Phi = Points(
        percolated_state_vectors,
        collect(1:size(percolated_state_vectors, 2)),
        fill(Inf, size(percolated_state_vectors, 2)),
        fill(Inf, size(percolated_state_vectors, 2)),
    )

    # ==========================================================================
    if verbose
        println("Calculating Survival Signature Entires...")
    end
    Phi.solution, Phi.confidence = Evaluation.computeSurvivalSignatureEntry(
        sys, sim, Phi.coordinates
    )
    if verbose
        println("Survival Signature Calculated\n")
    end
    # ==========================================================================

    # clean-up unnessesary method definitions
    methods = Methods(method, nothing, nothing, nothing, nothing, nothing, nothing)

    Phi = expandPhi!(Phi, state_vectors)        # resize Phi to full (non-percolated) version
    Phi = reshapePhi!(Phi)                      # reshape Phi to retangular arrays
    # for the purpose of comparison
    # might make more sense to start with this size
    # but many changes would be necessary   

    signature = Model(Phi, MonteCarloModel(), sys, sim, methods, nothing, nothing)

    # ==========================================================================
    # Timing
    elapsed_time = (time_ns() - start_time) / 1e9       # in seconds
    metrics = Metrics(elapsed_time)
    signature.metrics = metrics

    if verbose
        println("\tTime Elapsed: $(elapsed_time)s")
        println(" ")
        println("Finished Successfully.\n")
    end

    return signature
end

function simulate(
    method::IntervalPredictorSimulation,
    sys::System,
    sim::Simulation,
    methods::Methods;
    verbose::Bool=false,
    shape_parameter::Union{Nothing,Float64}=nothing,
)::Model
    if verbose
        printDetails(sys, sim, methods)
    end

    start_time = time_ns()

    # =============================== PERCOLATION ==============================
    state_vectors, percolated_state_vectors, sim.threshold = Evaluation.generateStateVectors(
        sys
    )

    Phi = Points(
        percolated_state_vectors,
        collect(1:size(percolated_state_vectors, 2)),
        fill(Inf, size(percolated_state_vectors, 2)),
        fill(Inf, size(percolated_state_vectors, 2)),
    )

    # ============================= STARTING POINTS ============================
    if verbose
        println("Generating Starting Points...")
    end
    starting_points = StartingPoints.generateStartingPoints(
        methods.starting_points_method, Phi.coordinates, sys.types
    )

    starting_points.solution, starting_points.confidence = Evaluation.computeSurvivalSignatureEntry(
        sys, sim, starting_points.coordinates
    )
    if verbose
        println("Starting Points Generated.\n")
    end

    # ================================ CENTERS =================================
    if verbose
        println("Generating Centers...")
    end
    centers = Centers.generateCenters(methods.centers_method, state_vectors, sim.threshold)
    if verbose
        println("Centers Generated.\n")
    end
    # ============================== CONSTRAINTS ===============================

    constraints = monotonicity_constraints(centers)

    # ============================ SHAPE PARAMETER =============================
    if verbose
        println("Compute Shape Parameters...")
    end
    if isnothing(shape_parameter)
        shape_parameter = ShapeParameter.computeShapeParameter(
            methods.shape_parameter_method, Phi.coordinates, starting_points, centers
        )
    else
        shape_parameter = shape_parameter
    end
    if verbose
        println("\tShape Parameter: $(shape_parameter)")
        println("Shape Parameter Computed.\n")
    end
    # ============================ BASIC FUNCTION ==============================
    if verbose
        println("Initializing Basis Function...")
    end
    starting_basis = BasisFunction.basis(
        methods.basis_function_method, shape_parameter, starting_points.coordinates, centers
    )
    if verbose
        println("Basis Function Initialized.\n")
    end

    # ============================ INITIAL WEIGHTS =============================
    if verbose
        println("Initializing Weights...")
    end
    # convex method for minimizing the least-squared error of the weights
    initial_weights = lsqr(starting_basis, starting_points.solution, constraints)
    if verbose
        println("Weights Initialized.\n")
    end
    # ========================== ADAPTIVE REFINEMENT ===========================
    if verbose
        println("Beginning Adaptive Refinement...")
    end
    # Taylor-Expansion Based Adaptive Design - Mo et al.
    evaluated_points, weights, upper_bound, lower_bound = AdaptiveRefinement.adaptiveRefinement(
        methods.adaptive_refinement_method,
        Phi,
        starting_points,
        sys,
        sim,
        methods,
        initial_weights,
        centers,
        constraints,
        shape_parameter,
    )
    if verbose
        println("Adaptive Refinement Completed\n")
    end

    # ========================== INTERVAL PREDICTOR ============================
    ipm = IntervalPredictorModel.intervalPredictor(
        evaluated_points,
        upper_bound,
        lower_bound,
        centers,
        shape_parameter,
        weights,
        methods,
    )

    Phi = mergePoints!(Phi, evaluated_points)   # fill evaluated values into Phi
    Phi = expandPhi!(Phi, state_vectors)        # resize Phi to full (non-percolated) version
    Phi = reshapePhi!(Phi)                      # reshape Phi to retangular arrays
    #                                           # for the purpose of comparison
    #                                           # might make more sense to start with this size
    #                                           # but many changes would be necessary   

    signature = Model(Phi, ipm, sys, sim, methods, nothing, nothing) # struct

    # ====================== EVALUATE REMAINING POINTS =========================
    if verbose
        println("Evaluating Remaining Points...")
    end
    signature = Evaluation.evaluate(signature)
    if verbose
        println("Remaining Points Evaluated.\n")
    end
    # ==========================================================================

    # Timing
    elapsed_time = (time_ns() - start_time) / 1e9       # in seconds
    metrics = Metrics(elapsed_time)
    signature.metrics = metrics

    if verbose
        println("Time Elapsed: $(elapsed_time)s")
        println(" ")
        println("Finished Successfully.\n")
    end

    return signature
end

# ==============================================================================

# ==============================================================================
# relocate this function
function mergePoints!(Phi::Points, evaluated::Points)::Points
    for (i, idx) in enumerate(evaluated.idx)
        Phi.solution[idx] = evaluated.solution[i]
        Phi.confidence[idx] = evaluated.confidence[i]
    end
    return Phi
end

function expandPhi!(Phi::Points, full_state_vector::Array)::Points
    solutions = zeros(size(full_state_vector, 2))
    confidence = zeros(size(full_state_vector, 2))

    for (i, state_vector) in enumerate(eachcol(full_state_vector))
        # find if that state_vector (vector) is in Phi.coordinates
        idx = findfirst(x -> x == state_vector, eachcol(Phi.coordinates))

        if idx !== nothing
            solutions[i] = Phi.solution[idx]
            confidence[i] = Phi.confidence[idx]
        else
            # set solutions[i] and confidence[i] to zero
            solutions[i] = 0.0
            confidence[i] = 0.0
        end
    end

    return Points(
        full_state_vector, collect(1:size(full_state_vector, 2)), solutions, confidence
    )
end

function reshapePhi!(Phi::Points)::Points
    # Create zero arrays with appropriate dimensions
    dimensions = Int.(Phi.coordinates[:, end])
    coordinates = Array{Vector{Int},length(dimensions)}(undef, dimensions...)
    index = zeros(Int, dimensions...)
    solution = zeros(Float64, dimensions...)
    confidence = zeros(Float64, dimensions...)

    for idx in CartesianIndices(coordinates)
        # Convert CartesianIndex to a Tuple and then to a Vector
        idx_vector = collect(Tuple(idx))

        # Store the vector at the current index
        coordinates[idx] = idx_vector

        # Find the corresponding index in Phi.coordinates
        id = findfirst(x -> x == idx_vector, eachcol(Phi.coordinates))

        # Ensure id is not nothing before assignment
        if id !== nothing
            index[idx] = Phi.idx[id]
            solution[idx] = Phi.solution[id]
            confidence[idx] = Phi.confidence[id]
        end
    end

    Phi = Points(coordinates, index, solution, confidence)

    return Phi
end
# ==============================================================================

end