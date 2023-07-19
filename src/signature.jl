"""
Calculate the exact survival signature of `system` in regards to the given `types` and structure function `φ`.

This evaluates all possible combinations of network states and is time consuming for large networks. An optional
`preprocessor` can be used to exclude parts of the survival signature beforhand.
"""
function survivalsignature(
    system::Any,
    types::Dict{Int64,Array{Int64,1}},
    φ::Function,
    preprocessor=nothing,
)
    Φ, _ = prepare_survival_signature(types)

    preprocessor !== nothing && (preprocessor(Φ, system))

    indices = findall(x -> isinf(x), Φ)
    entries = @showprogress pmap(indices) do index
        return exactentry(index, system, types, φ)
    end

    for (index, entry) in zip(indices, entries)
        Φ[index] = entry
    end

    return Φ
end

function survivalsignature(
    system::Any,
    types::Dict{Int64,Array{Int64,1}},
    φ::Function,
    samples::Int64,
    limit::Float64=0.001,
    preprocessor=nothing,
)
    Φ, components_per_type = prepare_survival_signature(types)
    preprocessor !== nothing && (preprocessor(Φ, system))

    cov = zeros(size(Φ))

    indices = findall(x -> isinf(x), Φ)
    entries = @showprogress pmap(indices) do index
        n_c = numberofcombinations(components_per_type, index)
        if n_c <= samples
            return exactentry(index, system, types, φ), 0
        else
            return approximateentry(index, system, types, φ, samples, limit)
        end
    end

    for (index, entry) in zip(indices, entries)
        s, c = entry
        Φ[index] = s
        cov[index] = c
    end

    return Φ, cov
end

include("rbf/basis.jl")
include("rbf/opt.jl")

struct IPMSurvivalSignature
    X::AbstractMatrix
    C::AbstractMatrix
    σ::Vector{<:Real}
    wmax::Vector{<:Real}
    wmin::Vector{<:Real}
end

function IPMSurvivalSignature(
    system::Any,
    types::Dict{Int,Vector{Int}},
    φ::Function,
    ci::Vector{Int};
    samples::Integer=10000,
    covtol::Real=1e-3,
    wtol::Real=1e-3
)
    components_per_type = ones(Int, length(types))
    for (type, components) in types
        components_per_type[type] += length(components)
    end

    fc = percolation(system)
    threshold = sum(components_per_type .- 1) * (1 - fc)

    Ω = mapreduce(t -> [t...]', vcat, Iterators.product([1:c for c in components_per_type]...))
    Ω = [Ω[i, :] for i in 1:size(Ω, 1) if sum(Ω[i, :] .- 1) >= threshold]
    Ω = float.(hcat(Ω...))

    lb = minimum(Ω; dims=2)
    ub = maximum(Ω; dims=2)

    C = hcat([x for x in eachcol(Ω)]...)
    tree = KDTree(C)

    ranges = [range(0.0, 1.0; length=l) for l in fill(2, length(types))]
    Xn = mapreduce(t -> [t...], hcat, Iterators.product(ranges...))
    Xn = (Xn .* (ub .- lb) .+ lb)

    idx, _ = nn(tree, Xn)
    idx = unique(idx) # in case two have the same nearest neighbor

    Xn = C[:, idx]
    # fn = @showprogress "Initial Points" map(eachcol(Xn)) do x
    fn = map(eachcol(Xn)) do x
        entry = CartesianIndex(Int.(x)...)
        if (numberofcombinations(components_per_type, entry)) <= samples
            return exactentry(entry, system, types, φ), 0
        else
            return approximateentry(entry, system, types, φ, samples, covtol)
        end
    end

    cn = getindex.(fn, 2)
    fn = getindex.(fn, 1)

    ranges = [range(l, u, c) for (l, u, c) in zip(lb, ub, ci)]
    centers = hcat([[c...] for c in Iterators.product(ranges...) if sum(c .- 1) > threshold]...)

    Q = (getindex.(ranges, 2) .- getindex.(ranges, 1)) ./ 2

    P = gaussian.(distance(Xn, centers, Q))
    Pc = gaussian.(distance(C, centers, Q))

    @show size(P)
    @show size(fn)
    @show size(centers)

    w = lsqr(P, fn, centers)

    Lmax = sqrt(sum((ub .- lb) .^ 2))

    stop = 3
    while stop < 2
        tree = KDTree(Xn)

        candidate_idx = Not(idx)
        candidates = @view C[:, candidate_idx]

        i, D = nn(tree, candidates)

        function s(x)
            return (gaussian.(distance(x, centers, Q))*w)[1]
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

        entry = CartesianIndex(Int.(C[:, c])...)
        if (numberofcombinations(components_per_type, entry)) <= samples
            ftead, ctead = exactentry(entry, system, types, φ), 0
        else
            ftead, ctead = approximateentry(entry, system, types, φ, samples, covtol)
        end

        push!(fn, ftead)
        push!(cn, ctead)
        push!(idx, c)

        P = basis(Xn, centers, Q)
        w_old = w
        w = solve(P, fn, centers)

        @show norm(w_old - w)
        stop = norm(w_old - w) < wtol ? stop + 1 : 0
    end

    f_u = min.(fn .+ (fn .* cn), 1.0)
    f_l = max.(fn .- (fn .* cn), 0.0)

    w_u, w_l = ipm(Xn, f_u, f_l, centers, Q)

    return IPMSurvivalSignature(Xn, centers, Q, w_u, w_l)
end

function exactentry(index::CartesianIndex, system, types, φ)
    combinations = []
    for (type, components) ∈ types
        push!(combinations, subsets(components, index[type] - 1))
    end

    results = map(
        x -> φ(system, Iterators.flatten(x) |> collect),
        Iterators.product(combinations...),
    )

    return count(results) / length(results)
end

function approximateentry(
    index::CartesianIndex,
    system,
    types::Dict{Int,Array{Int,1}},
    φ::Function,
    samples::Int,
    limit::Float64,
)
    working = 0
    Φ = 0
    cov = Inf

    for n = 1:samples
        # Generate a random network state for the current index
        x = []
        for (type, components) ∈ types
            shuffle!(components)
            append!(x, components[1:index[type]-1])
        end

        # Evaluate the structure function and update the current approximation of Φ
        φ(system, x) && (working += 1; Φ = working / n)
        Φ > 0 && (cov = √((Φ - Φ^2) / n) / Φ)

        # Break if the cov is below the defined threshold
        Φ != 1 && cov < limit && break

    end

    return Φ, cov
end

function numberofcombinations(components_per_type::Array{Int}, index::CartesianIndex)
    combinations = 1
    for (c, i) in zip(components_per_type, Tuple(index))
        combinations *= BigCombinatorics.Binomial(c - 1, i - 1)
    end
    return combinations
end

function prepare_survival_signature(types)
    components_per_type = ones(Int, length(types))
    for (type, components) in types
        components_per_type[type] += length(components)
    end

    Φ = ones(components_per_type...) * Inf

    return Φ, components_per_type
end
