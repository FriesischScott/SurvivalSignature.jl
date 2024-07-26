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

function floyd_warshall!(D)
    n = size(D, 1)

    @inbounds for k in 1:n, i in 1:n, j in 1:n
        sum = D[i, k] + D[k, j]
        if sum < D[i, j]
            D[i, j] = sum
        end
    end

    return D
end

function efficiency(D)
    n = size(D, 1)

    floyd_warshall!(D)
    D[D .== 0] .= Inf

    D = D .^ -1

    D[diagind(D)] .= 0

    return 1 / (n * (n - 1)) * sum(D)
end
