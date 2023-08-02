using LinearAlgebra

include("./berlin-metro-adjacency.jl")

degrees = vec(sum(A; dims=2))

types = Dict(1 => findall(degrees .<= 2.0), 2 => findall(degrees .> 2.0))

replace!(A, 0.0 => Inf) # ! Floy-Warshall
