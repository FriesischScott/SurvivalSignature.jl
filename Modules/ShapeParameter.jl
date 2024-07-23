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
using ..Structures: Hardy, Franke, Kuo, Rippa, BehrensdorfShape
using ..Structures: Gaussian

#include("SurvivalSignatureUtils.jl")
using ..SurvivalSignatureUtils

#include("BasisFunction.jl")
using ..BasisFunction

# ==============================================================================

export computeShapeParameter, selectShapeParameterMethod

# ==============================================================================

function selectShapeParameterMethod(
    method::String,
    coordinates::AbstractArray,
    centers::AbstractArray,
    starting_points::Points,
    confidence_interval::AbstractVector,
)
    method = lowercase(method)
    if method == "hardy"
        return Hardy(coordinates)
    elseif method == "franke"
        return Franke(coordinates)
    elseif method == "kuo"
        return Kuo(coordinates)
    elseif method == "rippa"
        return Rippa(starting_points, centers)
    elseif method == "behrensdorf"
        lb = minimum(coordinates; dims=2)
        ub = maximum(coordinates; dims=2)

        return BehrensdorfShape(confidence_interval, ub, lb)
    else
        error("Unrecognized Method: $method")
    end
end

# ============================= METHODS ========================================
function computeShapeParameter(Method::Hardy)

    # knn(k=2) returns the 2 closest points, since the 1. is itself 
    _, d = NearestNeighbors.knn(NearestNeighbors.KDTree(Method.points), Method.points, 2)

    d = sum(sum(d))

    return 1 / (0.815 * (d / size(Method.points, 2)))
end

function computeShapeParameter(Method::Franke)

    # leads to really slow adaptive refinement (change values are very high)
    points = SurvivalSignatureUtils.ensure_column_array!(Method.points)

    N::Int = length(points[:, 1])

    # D represents the diameter of the smallest circle which encompasses
    # all points, or the furthest distance between 2 points.
    D = largest_distance(points)

    return D / (sqrt(N) * 0.8)
end

function computeShapeParameter(Method::Kuo)
    # modified version of `franke()`
    # also really show adaptive refinement (possible do to formula flip)
    points = SurvivalSignatureUtils.ensure_column_array!(Method.points)
    N::Int = length(points[:, 1])

    D = largest_distance(points)
    return D / (nsqrt(N, 4) * 0.8)
end

function computeShapeParameter(Method::Rippa)
    cost_function =
        ϵ -> costFunction(
            Method.starting_points.solution,
            Method.starting_points.coordinates,
            Method.centers,
            ϵ,
        )
    result = optimize(cost_function, 0.1, 10.0)

    ϵ_opt = Optim.minimizer(result)
    # f_val = Optim.minimum(result)

    return ϵ_opt
end

function computeShapeParameter(method::BehrensdorfShape)
    lb = method.lower
    ub = method.upper
    ci = method.confidence_interval

    ranges = [range(l, u, c) for (l, u, c) in zip(lb, ub, ci)]

    return (getindex.(ranges, 2) .- getindex.(ranges, 1)) ./ 2
end

# ============================= UTILS ==========================================

function distance_matrix(points::AbstractArray)
    # Ensure points are in column format (rows as points)
    if size(points, 1) < size(points, 2)
        points = points'
    end
    return pairwise(Euclidean(), points; dims=2)
end

function largest_distance(points::AbstractArray)
    D = distance_matrix(points)
    return maximum(D)
end

function nsqrt(x::Number, n::Int)
    return x^(1 / n)
end

function costFunction(
    solutions::Union{AbstractVector,Number},    # true solutions of starting points
    coordinates::AbstractArray,                   # starting points
    centers::AbstractArray,
    shape_parameter::Number,
)
    N = length(solutions)

    A, _ = BasisFunction.basis(Gaussian(), shape_parameter, coordinates, centers)

    inv_A = pinv(A)  # Compute the pseudoinverse of the RBF matrix

    errors = zero(solutions)
    for i in 1:N
        A_i = A[:, i]
        inv_A_i = inv_A[i, :]  # Select the ith row of the pseudoinverse matrix
        y_pred_i = solutions[i] - (dot(A_i, inv_A_i) / inv_A_i[i])
        errors[i] = solutions[i] - y_pred_i
    end

    return LinearAlgebra.norm(errors)^2 / N
end

# =============================================================================

end