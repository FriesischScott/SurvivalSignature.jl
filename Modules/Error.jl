module Error

# this module to store all functions relating to the different errors between 
# true and estimated/appropriated Survival Signature values.

# ==============================================================================

using LinearAlgebra
using Statistics

using ..SurvivalSignatureUtils
using ..Structures: Model
using ..BasisFunction

# ==============================================================================
export calculateError
# ==============================================================================

function calculateError(
    methods::Vector{String}, models::Vector{Model}, compare::Union{Model,AbstractArray}
)
    errors = Vector{Dict}(undef, length(models))

    for (i, model) in enumerate(models)
        _, error = Error.calculateError(methods, model, compare)

        errors[i] = error
        model.errors = error
    end

    return models, errors
end

function calculateError(
    methods::Vector{String}, model::Model, compare::Union{Model,AbstractArray}
)
    errors = Dict{String,Union{String,Float64}}()

    compare_solution = isa(compare, Model) ? compare.Phi.solution : compare

    for method in methods
        errors[method] = calculateError(method, model.Phi.solution, compare_solution)
    end

    println("shape_parameter: $(model.str_method.shape_parameter_method)")
    println("samples: $(model.sim.samples)")
    SurvivalSignatureUtils._print(errors)

    return model, errors
end

function calculateError(
    methods::Vector{String}, arr1::AbstractArray, compare::Union{Model,AbstractArray}
)
    errors = Dict{String,Union{String,Float64}}()

    compare_solution = isa(compare, Model) ? compare.Phi.solution : compare

    for method in methods
        errors[method] = calculateError(method, arr1, compare_solution)
    end

    errors["name"] = "unknown"

    SurvivalSignatureUtils._print(errors)

    return errors
end

function calculateError(method::String, arr1::AbstractArray, arr2::AbstractArray)
    @assert prod(size(arr1)) == prod(size(arr2))

    func = selectErrorMethod(method)

    return func(arr1, arr2)
end

function selectErrorMethod(method::String)
    methods_dict = Dict("rmse" => rmse, "rae" => rae, "norm" => _norm, "nrsme" => nrmse)

    if haskey(methods_dict, lowercase(method))
        return methods_dict[lowercase(method)]
    else
        error("Unsupported Method: $method")
    end
end

# ==============================================================================

function rmse(arr1::AbstractArray, arr2::AbstractArray)
    # root mean sqauared error:

    # rmse = sqrt( Σ (y_i - ŷ_i)^2 / N )s
    # where:
    #   y_i  = arr1 
    #   ŷ_i = arr2
    #   N = total number of elements

    return sqrt(sum((arr1 .- arr2) .^ 2) / prod(size(arr1)))
end

function rae(true_values::AbstractArray, predicted_values::AbstractArray)
    # relative absolute error:

    # rae = Σ |y_i - ŷ_i| / Σ |y_i - ȳ|
    # where:
    #   y_i  = true_values
    #   ŷ_i = predicted_values
    #   ȳ = mean of true_values

    function absoluteSum(arr1::AbstractArray, arr2::Union{Number,AbstractArray})
        return sum(abs.(arr1 .- arr2))
    end

    return absoluteSum(true_values, predicted_values) /
           absoluteSum(true_values, Statistics.mean(true_values))
end

function _norm(weights::AbstractArray, old_weights::AbstractArray)
    # _norm to avoid overlap with third-party functions
    return LinearAlgebra.norm(weights - old_weights)
end

function nrmse()
    # to be filled out using the stopping critera from Mo et al.
    return nothing
end

# ==============================================================================

end
