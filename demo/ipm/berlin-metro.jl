using Distributed

using Distributions
using JLD2
using Plots

addprocs(12; exeflags="--project")

@load "demo/data/berlin-three-types-0-5-10000-MC.jld2"

@everywhere begin
    using SurvivalSignature

    include("../models/berlin-metro/berlin-metro-three-types.jl")

    φ = efficiency(nodes, A, 0.5)
end

N = 1000
covtol = 1e-3
wtol = 1e-3

ci = [15, 15, 10]

signature = survivalsignature(A, types, φ, ci; samples=N, covtol=covtol, wtol=wtol)

distributions = Dict(
    1 => Exponential(1 / 0.3), 2 => Exponential(1 / 0.4), 3 => Exponential(1 / 0.4)
)

t = range(0, 1, 100)

P = reliability(t, distributions, Φ)
Pl, Pu = reliability(t, distributions, signature)

p = plot(t, P; lw=1.5, label="MC Approximation")

plot!(p, t, Pl; color="black", style=:dash, lw=1.5, label="IPM Envelopes $N")
plot!(p, t, Pu; color="black", style=:dash, lw=1.5, label=:none)
