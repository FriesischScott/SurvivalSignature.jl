using Distributions
using JLD2 # extra dependency
using LinearAlgebra
using SurvivalSignature
using Plots
using Random
using Printf

n = 15
m = 15
adj = gridnetwork(n, m)
types = Dict(1 => collect(1:2:(n * m)), 2 => collect(2:2:(n * m)))
φ = s_t_connectivity([1:(n * m);], [1], [n * m])

covtol = 1e-3
wtol = 1e-3

ci = [15, 15]

function rmse(signature::IPMSurvivalSignature, Φ::Matrix)
    con = SurvivalSignature.monotonicity_constraints(signature.ipm.c)
    P = SurvivalSignature.basis(signature.X, signature.ipm.c, signature.ipm.σ)
    w = SurvivalSignature.lsqr(P, signature.f, con)
    threshold = sum(signature.k .- 1) * (1 - signature.fc)

    e = Float64[]

    for idx in CartesianIndices(Φ)
        if sum([Tuple(idx)...] .- 1) < threshold
            continue
        end
        s = (SurvivalSignature.basis([Tuple(idx)...], signature.ipm.c, signature.ipm.σ) * w)[1]
        push!(e, s - Φ[idx])
    end
    return sqrt(mean(e .^ 2))
end

function rae(signature::IPMSurvivalSignature, Φ::Matrix)
    con = SurvivalSignature.monotonicity_constraints(signature.ipm.c)
    P = SurvivalSignature.basis(signature.X, signature.ipm.c, signature.ipm.σ)
    w = SurvivalSignature.lsqr(P, signature.f, con)
    threshold = sum(signature.k .- 1) * (1 - signature.fc)

    e = Float64[]
    y = Float64[]

    for idx in CartesianIndices(Φ)
        if sum([Tuple(idx)...] .- 1) < threshold
            continue
        end
        s = (SurvivalSignature.basis([Tuple(idx)...], signature.ipm.c, signature.ipm.σ) * w)[1]
        push!(e, s - Φ[idx])
        push!(y, Φ[idx])
    end
    return mean(abs.(e)) / mean(abs.(y .- mean(y)))
end

signature_100 = IPMSurvivalSignature(
    adj, types, φ, ci; samples=10^2, covtol=covtol, wtol=wtol
)

signature_1000 = IPMSurvivalSignature(
    adj, types, φ, ci; samples=10^3, covtol=covtol, wtol=wtol
)

signature_10000 = IPMSurvivalSignature(
    adj, types, φ, ci; samples=10^4, covtol=covtol, wtol=wtol
)

@load "demo/data/grid-network-15x15-MC-10000.jld2"

println("+---------------------------------+")
println("| N     | M   | RMSE    |  RAE    |")
println("+---------------------------------+")
println(
    "| 100   | $(size(signature_100.X, 2)) | $(@sprintf("%.5f", rmse(signature_100, Φ))) | $(@sprintf("%.5f", rae(signature_100, Φ))) |",
)
println(
    "| 1000  | $(size(signature_1000.X, 2)) | $(@sprintf("%.5f", rmse(signature_1000, Φ))) | $(@sprintf("%.5f", rae(signature_1000, Φ))) |",
)
println(
    "| 10000 | $(size(signature_10000.X, 2)) | $(@sprintf("%.5f", rmse(signature_10000, Φ))) | $(@sprintf("%.5f", rae(signature_10000, Φ))) |",
)
println("+---------------------------------+")
