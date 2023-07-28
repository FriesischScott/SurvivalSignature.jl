function s_t_connectivity(nodes::Vector{Int}, source::Vector{Int}, target::Vector{Int})
    φ = function (system::AbstractMatrix, x::Vector)
        A = copy(system)
        s = copy(source)
        t = copy(target)

        # Remove failed nodes from s in case start and end nodes are the same
        s = filter(sᵢ -> sᵢ in x, s)

        failed = setdiff(nodes, x)

        A[failed, :] .= 0
        A[:, failed] .= 0

        while true
            if length(intersect(s, t)) > 0
                return true
            end

            A[:, s] .= 0
            s = findall(x -> x > 0, vec(sum(A[s, :]; dims=1)))

            if length(s) == 0
                return false
            end
        end
    end

    return φ
end

function efficiency(nodes::Vector{Int}, A::AbstractMatrix, loss::Real)
    E = efficiency(copy(A))

    return function φ(system::AbstractMatrix, x::Vector)
        A = copy(system)

        failed = setdiff(nodes, x)

        @inbounds A[failed, :] .= Inf
        @inbounds A[:, failed] .= Inf

        return efficiency(A) / E > loss
    end
end
