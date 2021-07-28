function s_t_connectivity(nodes::Vector{Int}, source::Vector{Int}, target::Vector{Int})

    φ = function (system::Array{Float64,2}, x::Vector)
        A = copy(system)
        s = copy(source)
        t = copy(target)

        failed = setdiff(nodes, x)

        A[failed, :] .= 0
        A[:, failed] .= 0

        while true
            if length(intersect(s, t)) > 0
                return true
            end

            A[:, s] .= 0
            s = findall(x -> x > 0, vec(sum(A[s, :], dims = 1)))

            if length(s) == 0
                return false
            end
        end

    end

    return φ
end
