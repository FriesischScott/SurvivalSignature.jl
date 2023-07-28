function gridnetwork(n::Int, m::Int)
    A = zeros(n * m, n * m)

    for i in 1:(n * m - 1)
        if !(i % n == 0)
            A[i, i + 1] = 1
        end
    end

    A[diagind(A, n)] .= 1

    return A += transpose(A)
end

function random_network(n::Int, p::Real)
    @assert n > 1
    @assert 0 <= p <= 1

    A = triu(rand(n, n) .<= p, 1)
    A += transpose(A)

    return float.(A)
end

function small_world_network(n::Int, k::Int, β::Real)
    @assert n > 1
    @assert iseven(k)

    A = zeros(Int, n, n)

    # connect to the k/2 neighbours on each side
    for i in 1:n, j in 1:n
        A[i, j] = 0 < abs((i - 1) - (j - 1)) % (n - 1 - k / 2) <= k / 2 ? 1 : 0
    end

    # rewire edge with probability β while avoiding self loops and duplicate edges
    for i in 1:n, j in 1:n
        if A[i, j] == 1 && rand() < β
            edges = findall(e -> e > 0, A[i, :])
            targets = setdiff(1:n, [i, edges...])
            A[i, j] = 0
            A[i, rand(targets)] = 1
        end
    end

    return float.(A)
end

function floyd_warshall!(D)
    replace!(D, 0 => Inf)
    n = size(D, 1)

    @inbounds for k in 1:n, i in 1:n, j in 1:n
        sum = D[i, k] + D[k, j]
        (sum < D[i, j]) && (D[i, j] = sum)
    end

    return D
end

function efficiency(D)
    n = size(D, 1)

    floyd_warshall!(D)
    D[D .== 0] .= Inf

    D = D .^ -1

    D[diagind(D)] == 0

    return 1 / (n * (n - 1)) * sum(D)
end
