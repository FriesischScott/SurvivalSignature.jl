function reliability(t::Vector{Float64}, Φ::Matrix{Float64}, dists::Dict)
    P = zeros(size(t))

    for index in CartesianIndices(Φ)
        p = ones(size(t))
        for k = 1:length(size(Φ))
            mk = size(Φ)[k] - 1
            lk = index[k] - 1
            Fk = dists[k]
            p .*=
                binomial(mk, lk) .* (cdf.(Fk, t) .^ (mk - lk)) .* ((1 .- cdf.(Fk, t)) .^ lk)
        end
        P .+= Φ[index] .* p
    end

    return P
end

function reliability(
    t::Vector{Float64},
    Φ::Matrix{Float64},
    types::Dict{Int64,Array{Int64,1}},
    failures::Matrix{Float64},
)
    idx = map_failures_to_type(failures, types)

    Φ_f = map_types_to_signature(idx, Φ, types)

    sort!(failures, dims = 2)

    P = zeros(size(t))

    for (r_i, r) ∈ enumerate(eachrow(failures))
        @inbounds P[t.<r[1]] .+= 1
        for c_i ∈ 1:length(r)-1
            if Φ_f[r_i, c_i] === 0.0
                break
            end
            @inbounds P[r[c_i].<=t.<r[c_i+1]] .+= Φ_f[r_i, c_i]
        end
        if Φ_f[r_i, end] > 0.0
            @inbounds P[t.>=r[end]] .+= Φ_f[r_i, end]
        end
    end

    return P ./ size(failures, 1)
end

function map_failures_to_type(f, types::Dict{Int64,Array{Int64,1}})
    type_map = Dict{Int64,Int64}()

    for n ∈ vcat(values(types)...)
        for (type, components) ∈ types
            n ∈ components && (type_map[n] = type; break)
        end
    end

    idx = zeros(size(f))

    for i ∈ 1:size(f, 1)
        @views idx[i, :] = map(x -> type_map[x], sortperm(f[i, :]))
    end

    return Int.(idx)
end

function map_types_to_signature(idx, Φ, types)
    F = ones(Int64, size(idx, 1), size(idx, 2), length(types))

    for (type, components) ∈ types
        F[:, :, type] += length(components) .- cumsum(idx .== type, dims = 2)
    end

    Φ_f = zeros(size(idx))

    for j = 1:size(F, 2)
        for i = 1:size(F, 1)
            Φ_f[i, j] = Φ[(F[i, j, :])...]
        end
    end

    return Φ_f
end
