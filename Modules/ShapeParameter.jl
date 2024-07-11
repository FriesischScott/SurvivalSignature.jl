module ShapeParameter

# ==============================================================================

using NearestNeighbors
using LinearAlgebra # needed for norm
using Optim
using Distances
using IterTools

# ==============================================================================

#include("Structures.jl")
using ..Structures: Points

#include("SurvivalSignatureUtils.jl")
using ..SurvivalSignatureUtils

#include("BasisFunction.jl")
using ..BasisFunction

# ============================= METHODS ========================================

function hardy(points::AbstractArray)

    # knn(k=2) returns the 2 closest points, since the 1. is itself 
    _, d = NearestNeighbors.knn(NearestNeighbors.KDTree(points), points, 2)

    d = sum(sum(d))

    return 1 / (0.815 * (d / size(points, 2)))
end

function franke(points::AbstractArray)

    # leads to really slow adaptive refinement (change values are very high)

    points = SurvivalSignatureUtils.ensure_column_array!(points)

    N::Int = length(points[:, 1])

    # D represents the diameter of the smallest circle which encompasses
    # all points, or the furthest distance between 2 points.
    D = largest_distance(points)

    return 0.8 * sqrt(N) / D
end

function kuo(points::Array)
    # modified version of `franke()`

    # also really show adaptive refinement (possible do to formula flip)

    N::Int = length(points[:, 1])

    D = largest_distance(points)

    return 0.8 * nsqrt(N, 4) / D
end

function rippa(X::AbstractArray, Y::Union{AbstractVector,Number}, centers::AbstractArray)
    # Utilizes LOOCV (leave-one-out cross-validation)

    cost_function = ϵ -> costFunction(Y, X, centers, ϵ)
    result = optimize(cost_function, 0.1, 10.0)

    ϵ_opt = Optim.minimizer(result)
    f_val = Optim.minimum(result)

    return ϵ_opt
end

function behrensdorf(ranges::AbstractArray)
    # spread in each direction (one per type)
    return (getindex.(ranges, 2) .- getindex.(ranges, 1)) ./ 2
end

# =========================== SELECTION ========================================

function selectShapeParameterMethod(method::String)
    methods_dict = Dict(
        "hardy" => hardy,
        "franke" => franke,
        "kuo" => kuo,
        "rippa" => rippa,
        "behrensdorf" => behrensdorf,
    )

    if haskey(methods_dict, lowercase(method))
        return methods_dict[lowercase(method)]
    else
        error("Unsupported Method: $method")
    end
end

function computeShapeParameter(
    method::String,
    points::AbstractArray,
    centers::AbstractArray,
    ranges::AbstractArray,
    partial_coordinates::AbstractArray,
    partial_solutions::Union{Nothing,AbstractArray,Number},
)
    func::Function = selectShapeParameterMethod(method)

    X::AbstractArray = partial_coordinates
    Y::Union{Nothing,AbstractArray,Number} = partial_solutions

    if func == rippa    # requires different inputs
        # X: starting points - Y: starting point outputs
        return rippa(X, Y, centers)
    elseif func == behrensdorf
        return behrensdorf(ranges)
    else
        return func(points)
    end
end

# ============================= UTILS ==========================================

function distance_matrix(points::AbstractArray)
    # Ensure points are in column format
    if size(points, 1) > size(points, 2)
        points = points'
    end
    return pairwise(Euclidean(), points; dims=2)
end

function largest_distance(points::AbstractArray)
    D = distance_matrix(points)
    return maximum(D)
end

function nearest_neighbor_distance(
    target_point::AbstractArray{<:Number}, points::AbstractArray
)
    if size(target_point, 1) < size(target_point, 2)
        target_point = target_point'
    end
    D = distance_matrix(hcat(target_point, points))
    return minimum(filter(x -> x > 0, D[:, 1]))
end

function nsqrt(x::Number, n::Int)
    return x^(1 / n)
end

function costFunction(
    y::Union{AbstractVector,Number},    # true solutions of starting points
    x::AbstractArray,                   # starting points
    centers::AbstractArray,
    shape_parameter::Number,
)
    N = length(y)

    A = BasisFunction.basis("gaussian", shape_parameter, x, centers)

    inv_A = pinv(A)  # Compute the pseudoinverse of the RBF matrix

    errors = zero(y)
    for i in 1:N
        A_i = A[:, i]
        inv_A_i = inv_A[i, :]  # Select the ith row of the pseudoinverse matrix
        y_pred_i = y[i] - (dot(A_i, inv_A_i) / inv_A_i[i])
        errors[i] = y[i] - y_pred_i
    end

    return LinearAlgebra.norm(errors)^2 / N
end

# =============================================================================

export *

end