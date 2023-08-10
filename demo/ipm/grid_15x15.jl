using Distributions
using JLD2 # extra dependency
using LinearAlgebra
using SurvivalSignature
using Plots
using Random

n = 15
m = 15
adj = gridnetwork(n, m)
types = Dict(1 => collect(1:2:(n * m)), 2 => collect(2:2:(n * m)))
φ = s_t_connectivity([1:(n * m);], [1], [n * m])

covtol = 1e-3
wtol = 1e-3

ci = [20, 20]

signature_100 = IPMSurvivalSignature(
    adj, types, φ, ci; samples=10^2, covtol=covtol, wtol=wtol
)
jldsave("grid-network-$(n)x$(m)-IPM-100.jld2"; signature_100)

signature_1000 = IPMSurvivalSignature(
    adj, types, φ, ci; samples=10^3, covtol=covtol, wtol=wtol
)
jldsave("grid-network-$(n)x$(m)-IPM-1000.jld2"; signature_1000)

signature_10000 = IPMSurvivalSignature(
    adj, types, φ, ci; samples=10^4, covtol=covtol, wtol=wtol
)
jldsave("grid-network-$(n)x$(m)-IPM-10000.jld2"; signature_10000)
