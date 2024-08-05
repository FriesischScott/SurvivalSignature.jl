module StructureCompilation

using ..Structures: Structures, Methods, Simulation
using ..Structures:
    SimulationType,
    StartingMethod,
    CentersMethod,
    ErrorType,
    ShapeParameterMethod,
    BasisFunctionMethod,
    AdaptiveRefinementMethod

# ==============================================================================

export compileMethods, compileSimulation

# ================================ METHODS =====================================
function compileMethods(
    sim_method::SimulationType,
    starting_points_method::StartingMethod,
    centers_method::CentersMethod,
    weight_change_method::ErrorType,
    shape_parameter_method::ShapeParameterMethod,
    basis_function_method::BasisFunctionMethod,
    adaptive_refinement_method::AdaptiveRefinementMethod,
)
    return Structures.Methods(
        sim_method,
        starting_points_method,
        centers_method,
        weight_change_method,
        shape_parameter_method,
        basis_function_method,
        adaptive_refinement_method,
    )
end

function compileMethods(
    sim_method::Vector{SimulationType},
    starting_points_method::StartingMethod,
    centers_method::CentersMethod,
    weight_change_method::ErrorType,
    shape_parameter_method::ShapeParameterMethod,
    basis_function_method::BasisFunctionMethod,
    adaptive_refinement_method::AdaptiveRefinementMethod,
)
    methods = Vector{Methods}(undef, length(sim_method))
    for (i, method) in enumerate(sim_method)
        methods[i] = Structures.Methods(
            method,
            starting_points_method,
            centers_method,
            weight_change_method,
            shape_parameter_method,
            basis_function_method,
            adaptive_refinement_method,
        )
    end

    return methods
end

function compileMethods(
    sim_method::SimulationType,
    starting_points_method::Vector{StartingMethod},
    centers_method::CentersMethod,
    weight_change_method::ErrorType,
    shape_parameter_method::ShapeParameterMethod,
    basis_function_method::BasisFunctionMethod,#
    adaptive_refinement_method::AdaptiveRefinementMethod,
)
    methods = Vector{Methods}(undef, length(starting_points_method))
    for (i, method) in enumerate(starting_points_method)
        methods[i] = Structures.Methods(
            sim_method,
            method,
            centers_method,
            weight_change_method,
            shape_parameter_method,
            basis_function_method,
            adaptive_refinement_method,
        )
    end

    return methods
end

function compileMethods(
    sim_method::SimulationType,
    starting_points_method::StartingMethod,
    centers_method::Vector{CentersMethod},
    weight_change_method::ErrorType,
    shape_parameter_method::ShapeParameterMethod,
    basis_function_method::BasisFunctionMethod,
    adaptive_refinement_method::AdaptiveRefinementMethod,
)
    methods = Vector{Methods}(undef, length(centers_method))
    for (i, method) in enumerate(centers_method)
        methods[i] = Structures.Methods(
            sim_method,
            starting_points_method,
            method,
            weight_change_method,
            shape_parameter_method,
            basis_function_method,
            adaptive_refinement_method,
        )
    end

    return methods
end

function compileMethods(
    sim_method::SimulationType,
    starting_points_method::StartingMethod,
    centers_method::Vector{CentersMethod},
    weight_change_method::ErrorType,
    shape_parameter_method::ShapeParameterMethod,
    basis_function_method::BasisFunctionMethod,
    adaptive_refinement_method::AdaptiveRefinementMethod,
)
    methods = Vector{Methods}(undef, length(weight_change_method))
    for (i, method) in enumerate(weight_change_method)
        methods[i] = Structures.Methods(
            sim_method,
            starting_points_method,
            centers_method,
            method,
            shape_parameter_method,
            basis_function_method,
            adaptive_refinement_method,
        )
    end

    return methods
end

function compileMethods(
    sim_method::SimulationType,
    starting_points_method::StartingMethod,
    centers_method::CentersMethod,
    weight_change_method::ErrorType,
    shape_parameter_method::Vector{ShapeParameterMethod},
    basis_function_method::BasisFunctionMethod,
    adaptive_refinement_method::AdaptiveRefinementMethod,
)
    methods = Vector{Methods}(undef, length(shape_parameter_method))
    for (i, method) in enumerate(shape_parameter_method)
        methods[i] = Structures.Methods(
            sim_method,
            starting_points_method,
            centers_method,
            weight_change_method,
            method,
            basis_function_method,
            adaptive_refinement_method,
        )
    end

    return methods
end

function compileMethods(
    sim_method::SimulationType,
    starting_points_method::StartingMethod,
    centers_method::CentersMethod,
    weight_change_method::ErrorType,
    shape_parameter_method::ShapeParameterMethod,
    basis_function_method::Vector{BasisFunctionMethod},
    adaptive_refinement_method::AdaptiveRefinementMethod,
)
    methods = Vector{Methods}(undef, length(basis_function_method))
    for (i, method) in enumerate(basis_function_method)
        methods[i] = Structures.Methods(
            sim_method,
            starting_points_method,
            centers_method,
            weight_change_method,
            shape_parameter_method,
            method,
            adaptive_refinement_method,
        )
    end

    return methods
end

function compileMethods(
    sim_method::SimulationType,
    starting_points_method::StartingMethod,
    centers_method::CentersMethod,
    weight_change_method::ErrorType,
    shape_parameter_method::ShapeParameterMethod,
    basis_function_method::BasisFunctionMethod,
    adaptive_refinement_method::Vector{AdaptiveRefinementMethod},
)
    methods = Vector{Methods}(undef, length(adaptive_refinement_method))
    for (i, method) in enumerate(adaptive_refinement_method)
        methods[i] = Structures.Methods(
            sim_method,
            starting_points_method,
            centers_method,
            weight_change_method,
            shape_parameter_method,
            basis_function_method,
            method,
        )
    end

    return methods
end

# ============================== SIMULATION ====================================

function compileSimulation(samples::Int, covtol::Float64, wtol::Float64)
    return Structures.Simulation(samples, covtol, wtol, nothing)
end

function compileSimulation(samples::Vector{Int}, covtol::Float64, wtol::Float64)
    sims = Vector{Simulation}(undef, length(samples))
    for (i, sample) in enumerate(samples)
        sims[i] = Structures.Simulation(sample, covtol, wtol, nothing)
    end

    return sims
end

function compileSimulation(samples::Int, covtols::Vector{Float64}, wtol::Float64)
    sims = Vector{Simulation}(undef, length(covtols))
    for (i, covtol) in enumerate(covtols)
        sims[i] = Structures.Simulation(samples, covtol, wtol, nothing)
    end

    return sims
end

function compileSimulation(samples::Int, covtol::Float64, wtols::Vector{Float64})
    sims = Vector{Simulation}(undef, length(wtols))
    for (i, wtol) in enumerate(wtols)
        sims[i] = Structures.Simulation(samples, covtol, wtol, nothing)
    end

    return sims
end

# ==============================================================================

end