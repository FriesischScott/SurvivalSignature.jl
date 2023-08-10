using Distributions
using JLD2 # extra dependency
using LinearAlgebra
using SurvivalSignature
using Plots
using Random

n = 15
m = 15
adj = gridnetwork(n, m)
types = Dict(1 => collect(3:2:((n * m) - 1)), 2 => collect(2:2:((n * m) - 1)))
φ = s_t_connectivity(
    [2:(n * m - 1);], findall(!iszero, adj[1, :]), findall(!iszero, adj[:, n * m])
)

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

# distributions = Dict(1 => Exponential(1), 2 => Weibull(2, 1))
# t = collect(range(0.3, 0.5, 200))

# Pl_1, Pu_1 = reliability(t, distributions, signature_1)
# Pl_2, Pu_2 = reliability(t, distributions, signature_2)
# Pl_3, Pu_3 = reliability(t, distributions, signature_3)

# P = reliability(t, distributions, Φ)
# p = plot(t, P; lw=1.5, label="MC Approximation")

# plot!(p, t, Pl_1; color="red", style=:dot, lw=1.5, label="IPM Envelopes 100")
# plot!(p, t, Pu_1; color="red", style=:dot, lw=1.5, label=:none)

# plot!(p, t, Pl_2; color="black", style=:dashdot, lw=1.5, label="IPM Envelopes 1000")
# plot!(p, t, Pu_2; color="black", style=:dashdot, lw=1.5, label=:none)

# plot!(p, t, Pl_3; color="green", style=:dash, lw=1.5, label="IPM Envelopes 10000")
# plot!(p, t, Pu_3; color="green", style=:dash, lw=1.5, label=:none)
