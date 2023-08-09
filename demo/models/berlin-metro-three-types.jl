using LinearAlgebra

include("./berlin-metro-adjacency.jl")

degrees = vec(sum(A; dims=2))

types = Dict(
    1 => findall(degrees .<= 2.0),
    2 => findall(2.0 .< degrees .<= 4.0),
    3 => findall(degrees .> 4.0),
)

replace!(A, 0.0 => Inf) # ! Floy-Warshall
A[diagind(A)] .= 0.0 # ! Floy-Warshall
