module Simulate

# this module to store the functions relating to simulating a given grid or model.
# ==============================================================================

using ..SurvivalSignatureUtils
using ..Structures: System, Simulation, StringMethods, Methods, Points, Model
using ..Structures: MonteCarloSimulation, RadialBasisSimulation, IntervalPredictorSimulation
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
    method::String,
    sys::System,
    sims::Vector{Simulation},
    str_methods::Vector{StringMethods},
)::Vector{Model}

    # Initialize the signatures array
    nsignatures = length(str_methods) * length(sims)
    signatures = Vector{Model}(undef, nsignatures)

    i = 1
    for str_method in str_methods
        for sim in sims
            signatures[i] = simulate(method, sys, sim, str_method)
            i += 1
        end
    end

    return signatures
end

function simulate(
    methods::Vector{String},
    sys::System,
    sims::Vector{Simulation},
    str_methods::Vector{StringMethods},
)::Vector{Model}

    # Length check
    len = length(methods)
    if length(str_methods) != len
        error("Simulation Parameters MUST be equal length.")
    end

    # Initialize the signatures array
    nsignatures = len * length(sims)
    signatures = Vector{Model}(undef, nsignatures)

    i = 1
    for (method, str_method) in zip(methods, str_methods)
        for sim in sims
            signatures[i] = simulate(method, sys, sim, str_method)
            i += 1
        end
    end

    return signatures
end

function simulate(
    methods::Vector{String},
    sys::System,
    sim::Simulation,
    str_methods::Vector{StringMethods},
)::Vector{Model}

    # Length check
    len = length(methods)
    if length(str_methods) != len
        error("Simulation Parameters MUST be equal length.")
    end

    # Initialize the signatures array
    signatures = Vector{Model}(undef, len)

    # Iterate over the elements of all vectors simultaneously
    for (i, (method, str_method)) in enumerate(zip(methods, str_methods))
        signatures[i] = simulate(method, sys, sim, str_method)
    end

    return signatures
end

function simulate(
    method::String, sys::System, sim::Simulation, str_methods::Vector{StringMethods}
)::Vector{Model}
    # Initialize the signatures array
    signatures = Vector{Model}(undef, length(str_methods))

    # Iterate over the elements of all vectors simultaneously
    for (i, str_method) in enumerate(str_methods)
        signatures[i] = simulate(method, sys, sim, str_method)
    end

    return signatures
end

function simulate(
    method::String, sys::System, sims::Vector{Simulation}, str_method::StringMethods
)::Vector{Model}

    # Initialize the signatures array
    signatures = Vector{Model}(undef, length(sims))

    # Iterate over the elements of all vectors simultaneously
    for (i, sim) in enumerate(sims)
        signatures[i] = simulate(method, sys, sim, str_method)
    end

    return signatures
end

function simulate(
    method::String, sys::System, sim::Simulation, str_method::StringMethods
)::Model
    # if the Method is given as a String - it is converted to the appropriate struct
    method = lowercase(method)

    struct_method = if method == "monte-carlo"
        MonteCarloSimulation()
    elseif method == "radial-basis-function"
        RadialBasisSimulation()
    elseif method == "interval-predictor"
        IntervalPredictorSimulation()
    else
        error("Unspecified Method: $method")
    end

    return simulate(struct_method, sys, sim, str_method)
end

function simulate(
    method::MonteCarloSimulation, sys::System, sim::Simulation, str_method::StringMethods
)::Model
    SurvivalSignatureUtils.header()
    printDetails(sys, sim, str_method)

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
    struct_method = Methods(method, nothing, nothing, nothing, nothing, nothing)

    Phi = expandPhi!(Phi, state_vectors)        # resize Phi to full (non-percolated) version
    Phi = reshapePhi!(Phi)                      # reshape Phi to retangular arrays
    # for the purpose of comparison
    # might make more sense to start with this size
    # but many changes would be necessary   

    signature = Model(Phi, nothing, sys, sim, struct_method, str_method, nothing)

    # ==========================================================================

    println("Finished Successfully.\n")

    return signature
end

function simulate(
    method::RadialBasisSimulation, sys::System, sim::Simulation, str_method::StringMethods
)::Model
    SurvivalSignatureUtils.header()
    printDetails(sys, sim, str_method)

    return signature
end

function simulate(
    method::IntervalPredictorSimulation,
    sys::System,
    sim::Simulation,
    str_method::StringMethods,
)::Model
    SurvivalSignatureUtils.header()
    printDetails(sys, sim, str_method)

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
    starting_points, starting_points_method_struct = StartingPoints.generateStartingPoints(
        str_method.starting_points_method, Phi.coordinates, sys.types
    )

    starting_points.solution, starting_points.confidence = Evaluation.computeSurvivalSignatureEntry(
        sys, sim, starting_points.coordinates
    )
    println("Starting Points Generated.\n")

    # ================================ CENTERS =================================
    println("Generating Centers...")
    centers, centers_method_struct = Centers.generateCenters(
        str_method.centers_method, state_vectors, sim.threshold, sim.confidence_interval
    )
    println("Centers Generated.\n")
    # ============================== CONSTRAINTS ===============================

    constraints = monotonicity_constraints(centers)

    # ============================ SHAPE PARAMETER =============================
    println("Compute Shape Parameters...")
    shape_parameter_struct = ShapeParameter.selectShapeParameterMethod(
        str_method.shape_parameter_method,
        Phi.coordinates,
        centers,
        starting_points,
        sim.confidence_interval,
    )

    shape_parameter = ShapeParameter.computeShapeParameter(shape_parameter_struct)
    println("\t$(str_method.shape_parameter_method): $shape_parameter")
    println("Shape Parameter Computed.\n")

    # ============================ BASIC FUNCTION ==============================
    println("Initializing Basis Function...")
    starting_basis, basis_function_struct = BasisFunction.basis(
        str_method, shape_parameter, starting_points.coordinates, centers
    )
    println("Basis Function Initialized.\n")

    # ============================ INITIAL WEIGHTS =============================
    println("Initializing Weights...")
    # convex method for minimizing the least-squared error of the weights
    initial_weights = lsqr(starting_basis, starting_points.solution, constraints)
    println("Weights Initialized.\n")
    # ========================== ADAPTIVE REFINEMENT ===========================

    method = Methods(
        method,
        starting_points_method_struct,
        centers_method_struct,
        str_method.weight_change_method,
        shape_parameter_struct,
        basis_function_struct,
    )

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

    signature = Model(Phi, ipm, sys, sim, method, str_method, nothing) # struct

    # ====================== EVALUATE REMAINING POINTS =========================
    println("Evaluating Remaining Points...")
    signature = Evaluation.evaluate(signature)
    println("Remaining Points Evaluated.\n")
    # ==========================================================================

    println("Finished Successfully.\n")

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