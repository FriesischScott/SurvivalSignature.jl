function percolation_preprocessor!(Φ, A)
    fc = percolation(A)

    threshold = sum(size(Φ) .- 1) * (1 - fc)

    for idx in CartesianIndices(Φ)
        if sum(Tuple(idx) .- 1) < threshold
            Φ[idx] = 0
        end
    end

    return nothing
end
