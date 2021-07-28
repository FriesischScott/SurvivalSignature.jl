using Distributed
using ClusterManagers
using JLD2 # extra dependency
using Distributions

addprocs(32, exeflags = "--project")

@everywhere using SurvivalSignature, LinearAlgebra

@everywhere joinpath(pwd(), "demo", "models", "reduced_gb.jl") |> include

@everywhere function floyd_warshall!(D)

    n = size(D, 1)

    @inbounds for k = 1:n, i = 1:n, j = 1:n
        sum = D[i, k] + D[k, j]
        (sum < D[i, j]) && (D[i, j] = sum)
    end

    return D
end

@everywhere function efficiency(D)
    n = size(D, 1)

    floyd_warshall!(D)
    D[D.==0] .= Inf

    D = D .^ -1

    D[diagind(D)] == 0

    1 / (n * (n - 1)) * sum(D)
end

@everywhere E = efficiency(copy(adj))

@everywhere function φ(system::Array{Float64,2}, x::Vector)

    A = copy(system)

    failed = setdiff(nodes, x)

    @inbounds A[failed, :] .= Inf
    @inbounds A[:, failed] .= Inf

    efficiency(A) / E > 0.5
end

@load "demo/data/reduced_gb_exact.jld2"

samples = [1e1, 1e2, 1e3, 1e4, 1e5, 1e6]
n = 1000
absolute_error = zeros(length(samples), n)
relative_error = zeros(length(samples), n)
mse = zeros(length(samples), n)

for (i, s) ∈ enumerate(samples)
    for j = 1:n
        signature, _ =
            survivalsignature(adj, types, φ, Int(s), 0.001, percolation_preprocessor)

        absolute_error[i, j] = norm(Φ - signature)
        relative_error[i, j] = norm(Φ - signature) / norm(Φ)
        mse[i, j] = mean((Φ - signature) .^ 2)
    end
end

@save "reduced_gb_errors.jld2" samples absolute_error relative_error mse
