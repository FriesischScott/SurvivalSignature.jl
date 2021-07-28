using Distributions
using JLD2 # extra dependency
using LinearAlgebra
using SurvivalSignature

adj = gridnetwork(6, 6)
types = Dict(1 => collect(1:2:36), 2 => collect(2:2:36))
φ = s_t_connectivity([1:36;], [1], [36])

# Load exact solution
@load "demo/data/grid_66_exact.jld2"

## Percolation error
signature_percolation = percolation_preprocessor(signature_66_exact, adj)

absolute_error_percolation = norm(signature_66_exact .- signature_percolation)
relative_error_percolation = absolute_error_percolation / norm(signature_66_exact)

# Approximation error
signature_monte_carlo, _ =
    survivalsignature(adj, types, φ, 10000, 0.0, percolation_preprocessor)

absolute_error_monte_carlo = norm(signature_percolation .- signature_monte_carlo) # Don't repeat the percolation error
relative_error_monte_carlo = absolute_error_monte_carlo / norm(signature_percolation)

# Reliability error

dists = Dict(1 => Exponential(1), 2 => Weibull(2, 1))
t = [0:0.001:1;]

P_exact = reliability(t, signature_66_exact, dists)
P_monte_carlo = reliability(t, signature_monte_carlo, dists)

absolute_reliability_error = norm(P_exact .- P_monte_carlo)
relative_reliability_error = absolute_reliability_error / norm(P_exact)


println("Absolute percolation error: $absolute_error_percolation")
println("Relative percolation error: $relative_error_percolation")

println("Absolute Monte Carlo error: $absolute_error_monte_carlo")
println("Relative Monte Carlo error: $relative_error_monte_carlo")

println("Absolute reliability error: $absolute_reliability_error")
println("Relative reliability error: $relative_reliability_error")
