module Error

# this module to store all functions relating to the different errors between 
# true and estimated/appropriated Survival Signature values.

# ==============================================================================

using LinearAlgebra
using Statistics

using ..SurvivalSignatureUtils

# ==============================================================================
export calculateError
# ==============================================================================

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

    # rmse = sqrt( Σ (y_i - ŷ_i)^2 / N )
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
