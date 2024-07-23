"""
Calculate the exact survival signature of `system` in regards to the given `types` and structure function `φ`.

This evaluates all possible combinations of network states and is time consuming for large networks. An optional
`preprocessor` can be used to exclude parts of the survival signature beforhand.
"""

using BigCombinatorics

function survivalsignature(
    system::Any, types::Dict{Int64,Array{Int64,1}}, φ::Function, preprocessor=nothing
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

function exactentry(index::CartesianIndex, system, types, φ)
    combinations = []
    for (type, components) in types
        push!(combinations, subsets(components, index[type] - 1))
    end

    results = map(
        x -> φ(system, collect(Iterators.flatten(x))), Iterators.product(combinations...)
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

    for n in 1:samples
        # Generate a random network state for the current index
        x = []
        for (type, components) in types
            shuffle!(components)
            append!(x, components[1:(index[type] - 1)])
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
