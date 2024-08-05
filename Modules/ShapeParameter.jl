module ShapeParameter

# ==============================================================================

using NearestNeighbors
using LinearAlgebra # needed for norm
using Optim
using Distances
using IterTools
using Plots
using Surrogates

# ==============================================================================

using ..Structures: Points
using ..Structures: Hardy, Rippa, DirectAMLS, IndirectAMLS
using ..Structures: Gaussian

using ..SurvivalSignatureUtils

using ..BasisFunction

# ==============================================================================

export computeShapeParameter

# ============================= METHODS ========================================
function computeShapeParameter(
    method::Hardy, points::Array, starting_points::Points, centers::Array
)::Float64

    # knn(k=2) returns the 2 closest points, since the 1. is itself 
    _, d = NearestNeighbors.knn(NearestNeighbors.KDTree(points), points, 2)

    d = sum(sum(d))

    return 1 / (0.815 * (d / size(points, 2)))
end

# 

function computeShapeParameter(
    method::Rippa, points::Array, starting_points::Points, centers::Array
)::Float64
    cost_function =
        ϵ -> optimized_LOOCV(
            starting_points.coordinates, starting_points.solution, centers, ϵ
        )
    result = optimize(cost_function, 0.1, 10.0)

    ϵ_opt = Optim.minimizer(result)
    # f_val = Optim.minimum(result)

    return ϵ_opt
end

function computeShapeParameter(
    method::DirectAMLS, points::Array, starting_points::Points, centers::Array
)
    cost_function =
        ϵ -> directAMLS(
            starting_points.coordinates,
            starting_points.solution,
            centers,
            ϵ,
            method.max_iterations,
            method.tolerance,
        )
    result = optimize(cost_function, 0.1, 10.0)

    ϵ_opt = Optim.minimizer(result)
    # f_val = Optim.minimum(result)

    return ϵ_opt
end

function computeShapeParameter(
    method::IndirectAMLS, points::Array, starting_points::Points, centers::Array
)
    cost_function =
        ϵ -> indirectAMLS(
            starting_points.coordinates,
            starting_points.solution,
            centers,
            ϵ,
            method.max_iterations,
            method.tolerance,
        )
    result = optimize(cost_function, 0.1, 10.0)

    ϵ_opt = Optim.minimizer(result)
    # f_val = Optim.minimum(result)

    return ϵ_opt
end

# ============================= UTILS ==========================================

function optimized_LOOCV(   # doesnt currently work
    coordinates::Array{Float64},
    solutions::Vector{Float64},
    centers::Array,
    shape_parameter::Float64,
)
    N = length(solutions)

    A = BasisFunction.basis(Gaussian(), shape_parameter, coordinates, centers)
    invA = pinv(A)

    errors = (invA * solutions) / diag(invA)

    return norm(errors)^2 / N
end

function directAMLS(   # doesnt currently work
    coordinates::Array{Float64},
    solutions::Vector{Float64},
    centers::Array,
    shape_parameter::Float64,
    max_iterations::Int,
    tolerance::Float64,
)

    # the Direct AMLS method necessitates A be a square matrix, however this is 
    # only the case if the number of centers matches the number of starting points.
    A = squareMatrix(BasisFunction.basis(Gaussian(), shape_parameter, coordinates, centers))
    I = Matrix(LinearAlgebra.I(size(A, 1)))

    v_prev = solutions # initializing
    M_prev = I
    cost_prev = Inf
    for n in 1:max_iterations
        v = (I - A) * v_prev .+ solutions

        # eigen dicomposition 
        eig = LinearAlgebra.eigen(I - A)
        Λ = eig.values
        X = eig.vectors

        M = Λ .* M_prev .+ I

        # cost vector
        e = v ./ diag(X * M * X')
        cost = LinearAlgebra.norm(e)

        # convergence
        if (cost - cost_prev) < tolerance
            return cost
        end

        cost_prev = cost
        M_prev = M
        v_prev = v
    end

    return cost
end

function squareMatrix(A::Matrix)::Matrix
    # the purpose of this function is to make a square-matrix from a non-square matrix
    if isSquare(A)
        # if the matrix is already square, this is unnessesary
        return A
    else
        return A * A'
    end
end

function isSquare(matrix::Matrix)::Bool
    return size(matrix, 1) == size(matrix, 2)
end

# continue reading Fasshauer and Zhang.
# there is a better method for this - Algorithm 6
function indirectAMLS(
    coordinates::Array{Float64},
    solutions::Vector{Float64},
    centers::Array,
    shape_parameter::Float64,
    max_iterations::Int,
    tolerance::Float64,
)
    N = length(solutions)

    error_prev = zeros(N)
    errors = zeros(N)

    for k in 1:N
        solution_k = solutions[k]
        coordinate_k = coordinates[:, k]

        Q = sum(
            solutions[j] * BasisFunction.basis(
                Gaussian(),
                shape_parameter,
                LinearAlgebra.norm(coordinate_k - coordinates[:, j]),
            ) for j in 1:N if j != k
        )

        for n in 1:max_iterations
            r = zeros(N)
            for j in 1:N
                r[j] = solutions[j] - Q
            end

            u = sum(
                r[j] * BasisFunction.basis(
                    Gaussian(),
                    shape_parameter,
                    LinearAlgebra.norm(coordinate_k - coordinates[:, j]),
                ) for j in 1:N if j != k
            )

            Q += u

            # error estimate
            errors[k] = abs(solution_k - Q)
        end

        cost = LinearAlgebra.norm(errors)

        if (cost - LinearAlgebra.norm(error_prev)) < tolerance
            return cost
        end
        error_prev .= errors
    end

    return cost
end

# ==============================================================================

end

# this function is slower, and therefore depreciated. however it is kept for comparison
# function computeShapeParameter(
#     method::Rippa2, points::Array, starting_points::Points, centers::Array
# )::Float64
#     cost_function =
#         ϵ -> LOOCV(
#             starting_points.coordinates, starting_points.solution, centers, ϵ
#         )
#     result = optimize(cost_function, 0.1, 5.0)

#     ϵ_opt = Optim.minimizer(result)
#     # f_val = Optim.minimum(result)

#     return ϵ_opt
# end
# function LOOCV(
#     coordinates::Array{Float64},
#     solutions::Vector{Float64},
#     centers::Array,
#     shape_parameter::Float64,
# )
#     N = length(solutions)

#     A = BasisFunction.basis(Gaussian(), shape_parameter, coordinates, centers)
#     invA = pinv(A)

#     errors = zeros(N)
#     for i in 1:N
#         Ai = @view A[:, i]
#         invAi = @view invA[i, :]
#         y_pred = solutions[i] - (dot(Ai, invAi) / invAi[i])
#         errors[i] = solutions[i] - y_pred
#     end
#     return norm(errors)^2 / N
# end