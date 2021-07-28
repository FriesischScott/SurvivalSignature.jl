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
