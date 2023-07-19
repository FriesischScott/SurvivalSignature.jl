using LinearAlgebra

include("./berlin-metro-adjacency.jl")

degrees = sum(A, dims=2) |> vec

types = Dict(
    1 => findall(degrees .<= 2.0),
    2 => findall(degrees .> 2.0),
)

# A[A.==0.0] .= Inf # ! Floy-Warshall

function floyd_warshall!(D)
    replace!(D, 0.0 => Inf)

    n = size(D, 1)

    @inbounds for k in 1:n, i in 1:n, j in 1:n
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

    D[diagind(D)] == 0

    1 / (n * (n - 1)) * sum(D)
end

E = efficiency(copy(A))

function φ(system::Array{Float64,2}, x::Vector)

    A = copy(system)

    failed = setdiff(nodes, x)

    @inbounds A[failed, :] .= Inf
    @inbounds A[:, failed] .= Inf

    return efficiency(A) / E > 0.5
end
