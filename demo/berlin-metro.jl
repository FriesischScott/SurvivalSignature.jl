using SurvivalSignature
using Distributions
using JLD2
using Plots

@load "demo/data/berlin-metro-signature-two-types-degree-2-0-5-10000.jld2"

# load adjency matrix, types and structure function
include("./models/berlin-metro-types.jl")

N = 1000
covtol = 1e-3
wtol = 1e-3

ci = [20, 20]

signature = IPMSurvivalSignature(A, types, φ, ci; samples=N, covtol=covtol, wtol=wtol)

distributions = Dict(1 => Exponential(1 / 0.25), 2 => Exponential(1 / 0.5))
t = [0.0:0.001:1.0;]

P = reliability(t, distributions, Φ)
Pl, Pu = reliability(t, distributions, signature)

p = plot(t, P; lw=1.5, label="Approximation")
plot!(p, t, Pl; color="black", style=:dash, lw=1.5, label="IPM Envelopes")
plot!(p, t, Pu; color="black", style=:dash, lw=1.5, label=:none)

# plot!(p, t, P4; color="red", style=:dash, lw=1.5, label="IPM Envelopes 2")
# plot!(p, t, P5; color="red", style=:dash, lw=1.5, label=:none)

# # plot!(p, t, P4; color="red", style=:dot, lw=1.5, label="MC Bounds")
# # plot!(p, t, P5; color="red", style=:dot, lw=1.5, label=:none)
