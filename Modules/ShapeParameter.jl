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
using ..Structures: Hardy, Franke, Kuo, Rippa
using ..Structures: Gaussian

#include("SurvivalSignatureUtils.jl")
using ..SurvivalSignatureUtils

#include("BasisFunction.jl")
using ..BasisFunction

# ==============================================================================

export computeShapeParameter

# ============================= METHODS ========================================
function computeShapeParameter(
    method::Hardy, points::Array, starting_points::Points, centers::Array
)

    # knn(k=2) returns the 2 closest points, since the 1. is itself 
    _, d = NearestNeighbors.knn(NearestNeighbors.KDTree(points), points, 2)

    d = sum(sum(d))

    return 1 / (0.815 * (d / size(points, 2)))
end

function computeShapeParameter(
    method::Franke, points::Array, starting_points::Points, centers::Array
)

    # leads to really slow adaptive refinement (change values are very high)
    points = SurvivalSignatureUtils.ensure_column_array!(points)

    N::Int = length(points[:, 1])

    # D represents the diameter of the smallest circle which encompasses
    # all points, or the furthest distance between 2 points.
    D = largest_distance(points)

    return D / (sqrt(N) * 0.8)
end

function computeShapeParameter(
    method::Kuo, points::Array, starting_points::Points, centers::Array
)
    # modified version of `franke()`
    # also really show adaptive refinement (possible do to formula flip)
    points = SurvivalSignatureUtils.ensure_column_array!(points)
    N::Int = length(points[:, 1])

    D = largest_distance(points)
    return D / (nsqrt(N, 4) * 0.8)
end

function computeShapeParameter(
    method::Rippa, points::Array, starting_points::Points, centers::Array
)
    cost_function =
        ϵ -> costFunction(starting_points.solution, starting_points.coordinates, centers, ϵ)
    result = optimize(cost_function, 0.1, 10.0)

    ϵ_opt = Optim.minimizer(result)
    # f_val = Optim.minimum(result)

    return ϵ_opt
end

# function computeShapeParameter(method::BehrensdorfShape)
#     lb = method.lower
#     ub = method.upper
#     ci = method.confidence_interval

#     ranges = [range(l, u, c) for (l, u, c) in zip(lb, ub, ci)]

#     return (getindex.(ranges, 2) .- getindex.(ranges, 1)) ./ 2
# end

# ============================= UTILS ==========================================

function distance_matrix(points::Array)
    # Ensure points are in column format (rows as points)
    if size(points, 1) < size(points, 2)
        points = points'
    end
    return pairwise(Euclidean(), points; dims=2)
end

function largest_distance(points::Array)
    D = distance_matrix(points)
    return maximum(D)
end

function nsqrt(x::Float64, n::Int)
    return x^(1 / n)
end

function costFunction(
    solutions::Union{Vector,Float64},    # true solutions of starting points
    coordinates::Array,                   # starting points
    centers::Array,
    shape_parameter::Float64,
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