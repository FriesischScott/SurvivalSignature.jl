function percolation_preprocessor(signature, adj)
    Φ = copy(signature)
    fc = percolation(adj)

    threshold = sum(size(Φ) .- 1) * (1 - fc)

    for idx in CartesianIndices(Φ)
        if sum(Tuple(idx) .- 1) < threshold
            Φ[idx] = 0
        end
    end

    return Φ
end
