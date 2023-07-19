using SurvivalSignature

# load adjency matrix, types and structure function
include("./models/berlin-metro-types.jl")

N = 1
covtol = 1e-3
wtol = 1e-3

ci = [20, 20]

IPMSurvivalSignature(A, types, φ, ci; samples=N, covtol=covtol, wtol=wtol)

# distributions = Dict(1 => Exponential(1 / 0.25), 2 => Exponential(1 / 0.5))
# t = [0.0:0.001:1.0;]

# function SurvivalSignature.reliability(t::Vector{Float64}, Φ::Matrix{Float64}, dists::Dict)
#     P = zeros(size(t))

#     for index in CartesianIndices(Φ)
#         p = ones(size(t))
#         for k = 1:length(size(Φ))
#             mk = size(Φ)[k] - 1
#             lk = index[k] - 1
#             Fk = dists[k]
#             p .*=
#                 BigCombinatorics.Binomial(mk, lk) .* (cdf.(Fk, t) .^ (mk - lk)) .* ((1 .- cdf.(Fk, t)) .^ lk)
#         end
#         P .+= Φ[index] .* p
#     end

#     return P
# end

# function reliability2(t::Vector{Float64}, Φ::Matrix{Float64}, dists::Dict, w)
#     P = zeros(size(t))

#     for index in CartesianIndices(Φ)
#         if sum([Tuple(index)...] .- 1) < threshold
#             continue
#         end
#         p = ones(size(t))
#         for k = 1:length(size(Φ))
#             mk = size(Φ)[k] - 1
#             lk = index[k] - 1
#             Fk = dists[k]
#             p .*=
#                 BigCombinatorics.Binomial(mk, lk) .* (cdf.(Fk, t) .^ (mk - lk)) .* ((1 .- cdf.(Fk, t)) .^ lk)
#         end
#         P .+= clamp(s([Tuple(index)...], w, centers, Q), 0.0, 1.0) .* p
#     end

#     return P
# end

# P = reliability(t, Φ, distributions)
# P2 = reliability2(t, Φ, distributions, w_l)
# P3 = reliability2(t, Φ, distributions, w_u)

# P4 = reliability2(t, Φ, distributions, w_l2)
# P5 = reliability2(t, Φ, distributions, w_u2)

# # P4 = reliability(t, Φ_u, distributions)
# # P5 = reliability(t, Φ_l, distributions)

# p = plot(t, P; lw=1.5, label="Approximation")
# plot!(p, t, P2; color="black", style=:dash, lw=1.5, label="IPM Envelopes")
# plot!(p, t, P3; color="black", style=:dash, lw=1.5, label=:none)

# plot!(p, t, P4; color="red", style=:dash, lw=1.5, label="IPM Envelopes 2")
# plot!(p, t, P5; color="red", style=:dash, lw=1.5, label=:none)

# # plot!(p, t, P4; color="red", style=:dot, lw=1.5, label="MC Bounds")
# # plot!(p, t, P5; color="red", style=:dot, lw=1.5, label=:none)
