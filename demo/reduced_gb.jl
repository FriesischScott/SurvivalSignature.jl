using SurvivalSignature
using LinearAlgebra
using JLD2 # extra dependency

joinpath(pwd(), "demo", "models", "reduced_gb.jl") |> include

function floyd_warshall!(D)

    n = size(D, 1)

    @inbounds for k = 1:n, i = 1:n, j = 1:n
        sum = D[i, k] + D[k, j]
        (sum < D[i, j]) && (D[i, j] = sum)
    end

    return D
end

function efficiency(D)
    n = size(D, 1)

    floyd_warshall!(D)
    D[D.==0] .= Inf

    D = D .^ -1

    D[diagind(D)] .= 0

    1 / (n * (n - 1)) * sum(D)
end

E = efficiency(copy(adj))

function φ(system::Array{Float64,2}, x::Vector)

    A = copy(system)

    failed = setdiff(nodes, x)

    @inbounds A[failed, :] .= Inf
    @inbounds A[:, failed] .= Inf

    efficiency(A) / E > 0.5
end

@load "demo/data/reduced_gb_exact.jld2"

signature, _ = survivalsignature(adj, types, φ, 10000, 0.001, percolation_preprocessor!)

absolute_error = norm(Φ - signature)
relative_error = norm(Φ - signature) / norm(Φ)

println("Absolute error: $absolute_error")
println("Relative error: $relative_error")
