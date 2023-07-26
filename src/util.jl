function gridnetwork(n::Int, m::Int)

    A = zeros(n * m, n * m)

    for i = 1:n*m-1
        if !(i % n == 0)
            A[i, i+1] = 1
        end
    end

    A[diagind(A, n)] .= 1

    A += transpose(A)
end

function random_network(n::Int, p::Real)
    @assert n > 1
    @assert 0 <= p <= 1

    A = triu(rand(n, n) .<= p, 1)
    return A + A'
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

    return A
end
