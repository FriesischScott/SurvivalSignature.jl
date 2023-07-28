using SurvivalSignature
using Distributions
using Plots
using LinearAlgebra
using Random
using JLD2

Random.seed!(39252)

n = 50
p = 0.2

A = random_network(n, p)

φ = efficiency([1:n;], A, 0.5)

nodes = rand(1:3, n)

types = Dict(
    1 => findall(nodes .== 1), 2 => findall(nodes .== 2), 3 => findall(nodes .== 3)
)

N = 1000
covtol = 1e-3
wtol = 1e-3

ci = [20, 20, 20]

distributions = Dict(
    1 => Exponential(1 / 0.25), 2 => Exponential(1 / 0.5), 3 => Exponential(1 / 0.25)
)
t = collect(range(0, 2, 1000))

signature = IPMSurvivalSignature(A, types, φ, ci; samples=N, covtol=covtol, wtol=wtol)
# Pl, Pu = reliability(t, distributions, signature)

Φ, cov = survivalsignature(A, types, φ, 1000)

jldsave("random-network.jld2"; signature, Φ, cov)

# P = reliability(t, distributions, Φ)
# p = plot(t, P; lw=1.5, label="MC Approximation")

# plot!(p, t, Pl; color="black", style=:dash, lw=1.5, label="IPM Envelopes $N")
# plot!(p, t, Pu; color="black", style=:dash, lw=1.5, label=:none)
