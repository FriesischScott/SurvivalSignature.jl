module Simulate

# this module to store the functions relating to simulating a given grid or model.
# ==============================================================================

using ..SurvivalSignatureUtils
using ..Structures: System, Simulation, Method, Points, Model
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
function selectSimulation(method::String)
    methods_dict = Dict(
        "monte-carlo" => monteCarloSimulation,
        "radial-basis-function" => radialBasisFunctionSimulation,
        "interval-predictor" => intervalPredictorSimulation,
    )

    if haskey(methods_dict, lowercase(method))
        return methods_dict[lowercase(method)]
    else
        error("Unsupported Method: $method")
    end
end

function simulate(sys::System, sim::Simulation, method::Method)
    func = selectSimulation(method.simulation_method)

    SurvivalSignatureUtils.header()
    printDetails(sys, sim, method)

    return func(sys, sim, method)
end
# ==============================================================================

# relocate this function
function mergePoints!(Phi::Points, evaluated::Points)
    for (i, idx) in enumerate(evaluated.idx)
        Phi.solution[idx] = evaluated.solution[i]
        Phi.confidence[idx] = evaluated.confidence[i]
    end
    return Phi
end

function expandPhi!(Phi::Points, full_state_vector::AbstractArray)
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

function reshapePhi!(Phi::Points)
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

function monteCarloSimulation(sys::System, sim::Simulation, method::Method)
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

    println("Calculating Survival Signature Entires...")

    Phi.solution, Phi.confidence = Evaluation.computeSurvivalSignatureEntry(
        sys, sim, Phi.coordinates
    )
    println("Survival Signature Calculated\n")

    # ==========================================================================

    # clean-up unnessesary method definitions
    method = Method("monte-carlo", "n/a", "n/a", "n/a", "n/a", "n/a", 1)

    Phi = expandPhi!(Phi, state_vectors)        # resize Phi to full (non-percolated) version
    Phi = reshapePhi!(Phi)                      # reshape Phi to retangular arrays
    # for the purpose of comparison
    # might make more sense to start with this size
    # but many changes would be necessary   

    signature = Model(Phi, nothing, sys, sim, method)

    # ==========================================================================

    println("Finished Successfully.\n")

    return signature
end

function radialBasisFunctionSimulation(sys::System, sim::Simulation, method::Method)
    return signature
end

function intervalPredictorSimulation(sys::System, sim::Simulation, method::Method)

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
    println("Generating Starting Points...")
    starting_points = StartingPoints.generateStartingPoints(
        method.starting_points_method, Phi.coordinates, sys.types
    )

    starting_points.solution, starting_points.confidence = Evaluation.computeSurvivalSignatureEntry(
        sys, sim, starting_points.coordinates
    )
    println("Starting Points Generated.\n")

    # ================================ CENTERS =================================
    println("Generating Centers...")
    centers, ranges = Centers.generateCenters(
        method.centers_method, Phi.coordinates, sim.threshold, sim.confidence_interval
    )
    println("Centers Generated.\n")
    # ============================== CONSTRAINTS ===============================

    constraints = monotonicity_constraints(centers)

    # ============================ SHAPE PARAMETER =============================
    println("Compute Shape Parameters...")
    shape_parameter = ShapeParameter.computeShapeParameter(
        method.shape_parameter_method,
        Phi.coordinates,
        centers,
        ranges,
        starting_points.coordinates,
        starting_points.solution,
    )
    println("Shape Parameters Computed.\n")

    # ============================ BASIC FUNCTION ==============================
    println("Initializing Basis Function...")
    starting_basis = BasisFunction.basis(
        method.basis_function_method,
        shape_parameter,
        starting_points.coordinates,
        centers,
        method.smoothness_factor,
    )
    println("Basis Function Initialized.\n")

    # ============================ INITIAL WEIGHTS =============================
    println("Initializing Weights...")
    # convex method for minimizing the least-squared error of the weights
    initial_weights = lsqr(starting_basis, starting_points.solution, constraints)
    println("Weights Initialized.\n")
    # ========================== ADAPTIVE REFINEMENT ===========================
    println("Beginning Adaptive Refinement...")

    # Taylor-Expansion Based Adaptive Design - Mo et al.
    evaluated_points, weights, upper_bound, lower_bound = AdaptiveRefinement.adaptiveRefinement(
        Phi,
        starting_points,
        sys,
        sim,
        method,
        initial_weights,
        centers,
        constraints,
        shape_parameter,
    )
    println("Adaptive Refinement Completed\n")

    # ========================== INTERVAL PREDICTOR ============================
    ipm = IntervalPredictorModel.intervalPredictor(
        evaluated_points,
        upper_bound,
        lower_bound,
        centers,
        shape_parameter,
        weights,
        method,
    )

    Phi = mergePoints!(Phi, evaluated_points)   # fill evaluated values into Phi
    Phi = expandPhi!(Phi, state_vectors)        # resize Phi to full (non-percolated) version
    Phi = reshapePhi!(Phi)                      # reshape Phi to retangular arrays
    # for the purpose of comparison
    # might make more sense to start with this size
    # but many changes would be necessary   

    signature = Model(Phi, ipm, sys, sim, method) # struct

    # ====================== EVALUATE REMAINING POINTS =========================
    println("Evaluating Remaining Points...")
    signature = Evaluation.evaluate(signature)
    println("Remaining Points Evaluated.\n")
    # ==========================================================================

    println("Finished Successfully.\n")

    return signature
end

# ==============================================================================

end