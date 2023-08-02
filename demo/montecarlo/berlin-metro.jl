using SurvivalSignature
using Distributions
using JLD2
using Plots

@load "demo/data/berlin-metro-signature-two-types-degree-2-0-5-10000.jld2"

# load adjency matrix, nodes and types
include("../models/berlin-metro-types.jl")

φ = efficiency(nodes, A, 0.5)

N = 100
covtol = 1e-3
wtol = 1e-3

ci = [20, 20]

distributions = Dict(1 => Exponential(1 / 0.25), 2 => Exponential(1 / 0.5))
t = [0.0:0.001:1.0;]

P = reliability(t, distributions, Φ)
p = plot(t, P; lw=1.5, label="MC Approximation")

signature = IPMSurvivalSignature(A, types, φ, ci; samples=N, covtol=covtol, wtol=wtol)
Pl, Pu = reliability(t, distributions, signature)

plot!(p, t, Pl; color="black", style=:dash, lw=1.5, label="IPM Envelopes $N")
plot!(p, t, Pu; color="black", style=:dash, lw=1.5, label=:none)
