using SurvivalSignature
using LinearAlgebra
using JLD2 # extra dependency

include(joinpath(pwd(), "demo", "models", "reduced_gb", "reduced_gb.jl"))

φ = efficiency(nodes, adj, 0.5)

@load "demo/data/reduced_gb_exact.jld2"

signature, _ = survivalsignature(adj, types, φ, 10000, 0.001, percolation_preprocessor!)

absolute_error = norm(Φ - signature)
relative_error = norm(Φ - signature) / norm(Φ)

println("Absolute error: $absolute_error")
println("Relative error: $relative_error")
