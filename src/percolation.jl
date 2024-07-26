function percolation(adjacency::AbstractMatrix)
    # Some structure functions need Inf instead of 0.0 for no connections
    A = copy(adjacency)
    replace!(A, Inf => 0.0)

    degrees = sum(A, dims=1)
    κ = mean(degrees .^ 2) / mean(degrees)

    return 1 - 1 / (κ - 1)
end
