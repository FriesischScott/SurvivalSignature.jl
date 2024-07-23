module StructureCompilation

using ..Structures: Structures, StringMethods, Simulation

# ==============================================================================

export compileMethods, compileSimulation

# ==============================================================================

# ================================ METHODS =====================================
function compileMethods(
    sim_method::String,
    starting_points_method::String,
    centers_method::String,
    weight_change_method::String,
    shape_parameter_method::String,
    basis_function_method::String,
    smoothness_factor::Int,
)
    if lowercase(shape_parameter_method) == "behrensdorf"
        basis_function_method = "behrensdorf"
    end

    return Structures.StringMethods(
        sim_method,
        starting_points_method,
        centers_method,
        weight_change_method,
        shape_parameter_method,
        basis_function_method,
        smoothness_factor,
    )
end

function compileMethods(
    sim_method::Vector{String},
    starting_points_method::String,
    centers_method::String,
    weight_change_method::String,
    shape_parameter_method::String,
    basis_function_method::String,
    smoothness_factor::Int,
)
    if lowercase(shape_parameter_method) == "behrensdorf"
        basis_function_method = "behrensdorf"
    end

    methods = Vector{StringMethods}(undef, length(sim_method))
    for (i, method) in enumerate(sim_method)
        methods[i] = Structures.StringMethods(
            method,
            starting_points_method,
            centers_method,
            weight_change_method,
            shape_parameter_method,
            basis_function_method,
            smoothness_factor,
        )
    end

    return methods
end

function compileMethods(
    sim_method::String,
    starting_points_method::Vector{String},
    centers_method::String,
    weight_change_method::String,
    shape_parameter_method::String,
    basis_function_method::String,
    smoothness_factor::Int,
)
    if lowercase(shape_parameter_method) == "behrensdorf"
        basis_function_method = "behrensdorf"
    end

    methods = Vector{StringMethods}(undef, length(starting_points_method))
    for (i, method) in enumerate(starting_points_method)
        methods[i] = Structures.StringMethods(
            sim_method,
            method,
            centers_method,
            weight_change_method,
            shape_parameter_method,
            basis_function_method,
            smoothness_factor,
        )
    end

    return methods
end

function compileMethods(
    sim_method::String,
    starting_points_method::String,
    centers_method::Vector{String},
    weight_change_method::String,
    shape_parameter_method::String,
    basis_function_method::String,
    smoothness_factor::Int,
)
    if lowercase(shape_parameter_method) == "behrensdorf"
        basis_function_method = "behrensdorf"
    end

    methods = Vector{StringMethods}(undef, length(centers_method))
    for (i, method) in enumerate(centers_method)
        methods[i] = Structures.StringMethods(
            sim_method,
            starting_points_method,
            method,
            weight_change_method,
            shape_parameter_method,
            basis_function_method,
            smoothness_factor,
        )
    end

    return methods
end

function compileMethods(
    sim_method::String,
    starting_points_method::String,
    centers_method::Vector{String},
    weight_change_method::String,
    shape_parameter_method::String,
    basis_function_method::String,
    smoothness_factor::Int,
)
    if lowercase(shape_parameter_method) == "behrensdorf"
        basis_function_method = "behrensdorf"
    end

    methods = Vector{StringMethods}(undef, length(weight_change_method))
    for (i, method) in enumerate(weight_change_method)
        methods[i] = Structures.StringMethods(
            sim_method,
            starting_points_method,
            centers_method,
            method,
            shape_parameter_method,
            basis_function_method,
            smoothness_factor,
        )
    end

    return methods
end

function compileMethods(
    sim_method::String,
    starting_points_method::String,
    centers_method::String,
    weight_change_method::String,
    shape_parameter_method::Vector{String},
    basis_function_method::String,
    smoothness_factor::Int,
)
    methods = Vector{StringMethods}(undef, length(shape_parameter_method))
    for (i, method) in enumerate(shape_parameter_method)
        if lowercase(method) == "behrensdorf"
            basis_function_method = "behrensdorf"
        end

        methods[i] = Structures.StringMethods(
            sim_method,
            starting_points_method,
            centers_method,
            weight_change_method,
            method,
            basis_function_method,
            smoothness_factor,
        )
    end

    return methods
end

function compileMethods(
    sim_method::String,
    starting_points_method::String,
    centers_method::String,
    weight_change_method::String,
    shape_parameter_method::String,
    basis_function_method::Vector{String},
    smoothness_factor::Int,
)
    methods = Vector{StringMethods}(undef, length(basis_function_method))
    for (i, method) in enumerate(basis_function_method)
        if lowercase(method) == "behrensdorf"
            method = "behrensdorf"
        end

        methods[i] = Structures.StringMethods(
            sim_method,
            starting_points_method,
            centers_method,
            weight_change_method,
            method,
            basis_function_method,
            smoothness_factor,
        )
    end

    return methods
end

# ============================== SIMULATION ====================================

function compileSimulation(samples::Int, covtol::Float64, wtol::Float64, ci::Vector{Int})
    return Structures.Simulation(samples, covtol, wtol, ci, nothing)
end

function compileSimulation(
    samples::Vector{Int}, covtol::Float64, wtol::Float64, ci::Vector{Int}
)
    sims = Vector{Simulation}(undef, length(samples))
    for (i, sample) in enumerate(samples)
        sims[i] = Structures.Simulation(sample, covtol, wtol, ci, nothing)
    end

    return sims
end

function compileSimulation(
    samples::Int, covtols::Vector{Float64}, wtol::Float64, ci::Vector{Int}
)
    sims = Vector{Simulation}(undef, length(covtols))
    for (i, covtol) in enumerate(covtols)
        sims[i] = Structures.Simulation(samples, covtol, wtol, ci, nothing)
    end

    return sims
end

function compileSimulation(
    samples::Int, covtol::Float64, wtols::Vector{Float64}, ci::Vector{Int}
)
    sims = Vector{Simulation}(undef, length(wtols))
    for (i, wtol) in enumerate(wtols)
        sims[i] = Structures.Simulation(samples, covtol, wtol, ci, nothing)
    end

    return sims
end

function compileSimulation(
    samples::Union{Int,Vector{Int}},
    covtol::Union{Float64,Vector{Float64}},
    wtol::Union{Float64,Vector{Float64}},
    ci::Vector{Int},
)
    println("this function was called.")

    return error()
end

# ==============================================================================

end