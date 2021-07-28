function percolation(adjacency::Array{Float64,2})

    degrees = sum(adjacency, dims = 1)
    κ = mean(degrees .^ 2) / mean(degrees)

    return 1 - 1 / (κ - 1)
end
