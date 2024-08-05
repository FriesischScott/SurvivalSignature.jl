module Error

# this module to store all functions relating to the different errors between 
# true and estimated/appropriated Survival Signature values.

# ==============================================================================

using LinearAlgebra
using Statistics

using ..SurvivalSignatureUtils
using ..Structures: Model
using ..Structures: ErrorType, RMSE, RAE, NORM, NRMSE
using ..BasisFunction

# ==============================================================================
export calculateError
# ==============================================================================

function calculateError(
    methods::Union{ErrorType,Vector{ErrorType}},
    models::Array{Model},
    compare::Union{Model,Array};
    verbose::Bool=false,
)
    errors = Array{Dict{String,Float64}}(undef, size(models))

    for (i, model) in enumerate(models)
        _, error = Error.calculateError(methods, model, compare; verbose=verbose)

        errors[i] = error
        model.errors = error
    end

    return models, errors
end

function calculateError(
    methods::Union{ErrorType,Vector{ErrorType}},
    model::Model,
    compare::Union{Model,Array};
    verbose::Bool=false,
)
    error = Dict{String,Float64}()

    compare_solution = isa(compare, Model) ? compare.Phi.solution : compare
    methods = isa(methods, Vector) ? methods : [methods]

    for method in methods
        error[string(nameof(typeof(method)))] = calculateError(
            method, model.Phi.solution, compare_solution
        )
    end

    if verbose
        println(
            "shape_parameter: $(string(nameof(typeof((model.method.shape_parameter_method)))))",
        )
        println("samples: $(model.sim.samples)")
        SurvivalSignatureUtils._print(error)
    end

    return model, error
end

function calculateError(
    methods::Union{ErrorType,Vector{ErrorType}},
    arr1::Array{Float64},           # difference to include arrays, rather than entire models
    compare::Union{Model,Array};
    verbose::Bool=false,
)
    errors = Dict{String,Float64}()

    compare_solution = isa(compare, Model) ? compare.Phi.solution : compare
    methods = isa(methods, Vector) ? methods : [methods]

    for method in methods
        errors[string(nameof(typeof(method)))] = calculateError(
            method, arr1, compare_solution
        )
    end

    errors["name"] = "unknown"

    if verbose
        SurvivalSignatureUtils._print(errors)
    end

    return errors
end

# ==============================================================================

function calculateError(method::RMSE, arr1::Array{Float64}, arr2::Array{Float64})
    # root mean sqauared error:

    # rmse = sqrt( Σ (y_i - ŷ_i)^2 / N )s
    # where:
    #   y_i  = arr1 
    #   ŷ_i = arr2
    #   N = total number of elements

    return sqrt(sum((arr1 .- arr2) .^ 2) / prod(size(arr1)))
end

function calculateError(
    method::RAE, true_values::Array{Float64}, predicted_values::Array{Float64}
)
    # relative absolute error:

    # rae = Σ |y_i - ŷ_i| / Σ |y_i - ȳ|
    # where:
    #   y_i  = true_values
    #   ŷ_i = predicted_values
    #   ȳ = mean of true_values

    function absoluteSum(arr1::Array{Float64}, arr2::Union{Number,Array{Float64}})
        return sum(abs.(arr1 .- arr2))
    end

    return absoluteSum(true_values, predicted_values) /
           absoluteSum(true_values, Statistics.mean(true_values))
end

function calculateError(method::NORM, weights::Array{Float64}, old_weights::Array{Float64})
    # _norm to avoid overlap with third-party functions
    return LinearAlgebra.norm(weights - old_weights)
end

function calculateError(method::NRMSE)
    # to be filled out using the stopping critera from Mo et al.
    return nothing
end

# ==============================================================================

end
