using JuMP
using JLD2
using LinearAlgebra
using StatsBase
using Plots
using Distributions
using BigCombinatorics
using SurvivalSignature
using Random
using SCS
using NearestNeighbors
using InvertedIndices
using ForwardDiff
using ProgressMeter

@load "berlin-metro-signature-two-types-degree-2-0-5-10000.jld2"

include("two-types.jl")

components_per_type = ones(Int, length(types))
for (type, components) in types
    components_per_type[type] += length(components)
end

plotlyjs()

# Prepare the inputs
A[A.==Inf] .= 0.0
fc = SurvivalSignature.percolation(A)
threshold = sum(size(Φ) .- 1) * (1 - fc)
A[A.==0.0] .= Inf

dim = [size(Φ)...]

Ω = mapreduce(t -> [t...]', vcat, Iterators.product(1:size(Φ, 1), 1:size(Φ, 2)))
Ω = [Ω[i, :] for i in 1:size(Ω, 1) if sum(Ω[i, :] .- 1) >= threshold]
Ω = float.(hcat(Ω...))

cov[cov.==Inf] .= 0.0

f = [Φ[Int.(x)...] for x in eachcol(Ω)]
c = [cov[Int.(x)...] for x in eachcol(Ω)]
c[c.==Inf] .= 0.0

function basis(X, Y, Q)
    P = [exp(-sum((x .- c) .^ 2 .* Q)) for (x, c) in Iterators.product(eachcol(X), eachcol(Y))]
    return P ./ sum(P; dims=2)
end

function solve(P, f, C)
    model = Model(SCS.Optimizer)
    set_silent(model)
    @variable(model, x[1:size(C, 2)])
    @variable(model, res[1:length(f)])
    @constraint(model, res .== P * x - f)

    for (i, ci) in enumerate(eachcol(C))
        for (j, cj) in enumerate(eachcol(C))
            if ci[1] < cj[1] && ci[2] == cj[2]
                @constraint(model, x[i] <= x[j])
            elseif ci[2] < cj[2] && ci[1] == cj[1]
                @constraint(model, x[i] <= x[j])
            end
        end
    end

    @objective(model, Min, sum(res .^ 2))

    JuMP.optimize!(model)
    return value.(x)
end

function s(x, w, C, Q)
    if sum(x .- 1) < threshold
        return 0.0
    end
    P = basis(x, C, Q)
    return sum(P * w)
end

function tead(Ω::AbstractMatrix, f::AbstractVector, cov::AbstractVector)
    lb = minimum(Ω; dims=2)
    ub = maximum(Ω; dims=2)

    C = hcat([x for x in eachcol(Ω)]...)
    tree = KDTree(C)

    ranges = [range(0.0, 1.0; length=l) for l in [5, 5]]
    Xn = mapreduce(t -> [t...], hcat, Iterators.product(ranges...))
    Xn = (Xn .* (ub .- lb) .+ lb)

    idx, _ = nn(tree, Xn)

    Xn = C[:, idx]
    fn = f[idx]
    # fn = @showprogress "Initial Points" map(eachcol(Xn)) do x
    #     entry = CartesianIndex(Int.(x)...)
    #     if (SurvivalSignature.numberofcombinations(components_per_type, entry)) <= 10_000
    #         return SurvivalSignature.exactentry(entry, A, types, φ)
    #     else
    #         return SurvivalSignature.approximateentry(entry, A, types, φ, 10_000, 1e-3)
    #     end
    # end

    Nx = 20
    Ny = 20

    x = range(lb[1], ub[1], Nx)
    y = range(lb[2], ub[2], Ny)

    centers = hcat([[c...] for c in Iterators.product(x, y) if sum(c .- 1) > threshold]...)

    σ1 = (x[2] - x[1]) ./ 2
    σ2 = (y[2] - y[1]) ./ 2

    Q = [1 / 2σ1^2, 1 / 2σ2^2]

    P = basis(Xn, centers, Q)
    Pc = basis(C, centers, Q)

    w = solve(P, fn, centers)

    Lmax = sqrt(sum((ub .- lb) .^ 2))
    # p = scatter(Xn[1, :], Xn[2, :], label="Initial")
    # display(p)

    stop = 0
    while stop < 2
        tree = KDTree(Xn)

        candidate_idx = Not(idx)
        candidates = @view C[:, candidate_idx]

        i, D = nn(tree, candidates)

        function s(x)
            return (basis(x, centers, Q)*w)[1]
        end

        ∇s = [zeros(size(Xn, 1)) for _ in 1:size(Xn, 2)]
        ∇s = map(ForwardDiff.DiffResults.GradientResult, ∇s)

        ∇s = map((r, x) -> ForwardDiff.gradient!(r, s, x), ∇s, eachcol(Xn))
        vs = ForwardDiff.DiffResults.value.(∇s)
        ∇s = ForwardDiff.DiffResults.gradient.(∇s)

        t = @views vs[i] + map((∇, a) -> dot(∇, a), ∇s[i], eachcol(candidates .- Xn[:, i]))

        R = @views abs.(Pc[candidate_idx, :] * w - t)

        J = D ./ maximum(D) + (1 .- D ./ Lmax) .* (R ./ maximum(R))

        _, c = findmax(J)
        Cn = candidates[:, c]

        c = findfirst(x -> all(x .== Cn), eachcol(Ω))

        Xn = hcat(Xn, C[:, c])
        push!(fn, f[c])
        push!(idx, c)

        # scatter!(p, [Xn[1, end]], [Xn[2, end]], label=:none)
        # display(p)

        P = basis(Xn, centers, Q)
        w_old = w
        w = solve(P, fn, centers)

        @show norm(w_old - w)
        stop = norm(w_old - w) < 1e-3 ? stop + 1 : 0
    end

    return Xn, fn, cov[idx], centers, w, Q
end

Xn, fn, cn, centers, w, Q = tead(Ω, f, c)

# fn = @showprogress "Initial Points" map(eachcol(Xn)) do x
#     entry = CartesianIndex(Int.(x)...)
#     if (SurvivalSignature.numberofcombinations(components_per_type, entry)) <= 1_000
#         return SurvivalSignature.exactentry(entry, A, types, φ)
#     else
#         return SurvivalSignature.approximateentry(entry, A, types, φ, 1_000, 1e-3)
#     end
# end

λ = 1.0

f_u = min.(fn .+ λ .* (fn .* cn), 1.0)
f_l = max.(fn .- λ .* (fn .* cn), 0.0)

function ipm(X::AbstractMatrix, f_u::AbstractVector, f_l::AbstractVector, centers::AbstractMatrix, Q::AbstractVector, w::AbstractVector)
    P = basis(X, centers, Q)

    nc = size(centers, 2)

    model = Model(SCS.Optimizer)
    @variable(model, x[1:nc*2])
    @variable(model, spread[1:length(f_u)])

    @constraint(model, spread .== P * (x[1:nc] - x[nc+1:end]))
    @objective(model, Min, mean(spread))

    @constraint(model, P * x[1:nc] .>= f_u)
    @constraint(model, P * x[nc+1:end] .<= f_l)

    # @constraint(model, x[1:nc] .>= w)
    # @constraint(model, w .>= x[nc+1:end])

    @constraint(model, x[1:nc] .>= x[nc+1:end])

    set_start_value.(x[1:nc], w)
    set_start_value.(x[nc+1], w)

    for (i, ci) in enumerate(eachcol(centers))
        for (j, cj) in enumerate(eachcol(centers))
            if ci[1] < cj[1] && ci[2] == cj[2]
                @constraint(model, x[i] <= x[j])
                @constraint(model, x[i+nc] <= x[j+nc])
            elseif ci[2] < cj[2] && ci[1] == cj[1]
                @constraint(model, x[i] <= x[j])
                @constraint(model, x[i+nc] <= x[j+nc])
            end
        end
    end

    JuMP.optimize!(model)

    w_u = value.(x[1:nc])
    w_l = value.(x[nc+1:end])
    return w_u, w_l
end

function ipm2(X::AbstractMatrix, Ω::AbstractMatrix, f_u::AbstractVector, f_l::AbstractVector, centers::AbstractMatrix, Q::AbstractVector, w::AbstractVector)
    P = basis(Ω, centers, Q)

    P2 = basis(X, centers, Q)

    nc = size(centers, 2)

    model = Model(SCS.Optimizer)
    @variable(model, x[1:nc*2])
    @variable(model, spread[1:size(Ω, 2)])

    @constraint(model, spread .== P * (x[1:nc] - x[nc+1:end]))
    @objective(model, Min, mean(spread))

    @constraint(model, P2 * x[1:nc] .>= f_u)
    @constraint(model, P2 * x[nc+1:end] .<= f_l)

    # @constraint(model, x[1:nc] .>= w)
    # @constraint(model, w .>= x[nc+1:end])

    @constraint(model, x[1:nc] .>= x[nc+1:end])

    set_start_value.(x[1:nc], w)
    set_start_value.(x[nc+1], w)

    for (i, ci) in enumerate(eachcol(centers))
        for (j, cj) in enumerate(eachcol(centers))
            if ci[1] < cj[1] && ci[2] == cj[2]
                @constraint(model, x[i] <= x[j])
                @constraint(model, x[i+nc] <= x[j+nc])
            elseif ci[2] < cj[2] && ci[1] == cj[1]
                @constraint(model, x[i] <= x[j])
                @constraint(model, x[i+nc] <= x[j+nc])
            end
        end
    end

    JuMP.optimize!(model)

    w_u = value.(x[1:nc])
    w_l = value.(x[nc+1:end])
    return w_u, w_l
end

w_u, w_l = ipm(Xn, f_u, f_l, centers, Q, w)

w_u2, w_l2 = ipm2(Xn, Ω, f_u, f_l, centers, Q, w)
# lb = minimum(Xn; dims=2)
# ub = maximum(Xn; dims=2)

# plot(lb[1]:ub[1], lb[2]:ub[2], (x, y) -> s([x; y], w, centers, Q); st=:surface)

distributions = Dict(1 => Exponential(1 / 0.25), 2 => Exponential(1 / 0.5))
t = [0.0:0.001:1.0;]

function SurvivalSignature.reliability(t::Vector{Float64}, Φ::Matrix{Float64}, dists::Dict)
    P = zeros(size(t))

    for index in CartesianIndices(Φ)
        p = ones(size(t))
        for k = 1:length(size(Φ))
            mk = size(Φ)[k] - 1
            lk = index[k] - 1
            Fk = dists[k]
            p .*=
                BigCombinatorics.Binomial(mk, lk) .* (cdf.(Fk, t) .^ (mk - lk)) .* ((1 .- cdf.(Fk, t)) .^ lk)
        end
        P .+= Φ[index] .* p
    end

    return P
end

function reliability2(t::Vector{Float64}, Φ::Matrix{Float64}, dists::Dict, w)
    P = zeros(size(t))

    for index in CartesianIndices(Φ)
        if sum([Tuple(index)...] .- 1) < threshold
            continue
        end
        p = ones(size(t))
        for k = 1:length(size(Φ))
            mk = size(Φ)[k] - 1
            lk = index[k] - 1
            Fk = dists[k]
            p .*=
                BigCombinatorics.Binomial(mk, lk) .* (cdf.(Fk, t) .^ (mk - lk)) .* ((1 .- cdf.(Fk, t)) .^ lk)
        end
        P .+= clamp(s([Tuple(index)...], w, centers, Q), 0.0, 1.0) .* p
    end

    return P
end

P = reliability(t, Φ, distributions)
P2 = reliability2(t, Φ, distributions, w_l)
P3 = reliability2(t, Φ, distributions, w_u)

P4 = reliability2(t, Φ, distributions, w_l2)
P5 = reliability2(t, Φ, distributions, w_u2)

# P4 = reliability(t, Φ_u, distributions)
# P5 = reliability(t, Φ_l, distributions)

p = plot(t, P; lw=1.5, label="Approximation")
plot!(p, t, P2; color="black", style=:dash, lw=1.5, label="IPM Envelopes")
plot!(p, t, P3; color="black", style=:dash, lw=1.5, label=:none)

plot!(p, t, P4; color="red", style=:dash, lw=1.5, label="IPM Envelopes 2")
plot!(p, t, P5; color="red", style=:dash, lw=1.5, label=:none)

# plot!(p, t, P4; color="red", style=:dot, lw=1.5, label="MC Bounds")
# plot!(p, t, P5; color="red", style=:dot, lw=1.5, label=:none)
