function lsqr(P::AbstractMatrix, f::AbstractVector, C::AbstractMatrix)
    model = Model(SCS.Optimizer)
    set_silent(model)
    @variable(model, x[1:size(C, 2)])
    @variable(model, res[1:length(f)])
    @constraint(model, res .== P * x - f)

    for (i, ci) in enumerate(eachcol(C)), (j, cj) in enumerate(eachcol(C))
        if sum(ci .< cj) == 1 && sum(ci .== cj) == (size(C, 1) - 1)
            @constraint(model, x[i] <= x[j])
        end
    end

    @objective(model, Min, sum(res .^ 2))

    JuMP.optimize!(model)
    return value.(x)
end

function ipm(X::AbstractMatrix, f_u::AbstractVector, f_l::AbstractVector, centers::AbstractMatrix, Q::AbstractVector)
    P = basis(X, centers, Q)

    nc = size(centers, 2)

    model = Model(SCS.Optimizer)
    @variable(model, x[1:nc*2])
    @variable(model, spread[1:length(f_u)])

    @constraint(model, spread .== P * (x[1:nc] - x[nc+1:end]))
    @objective(model, Min, mean(spread))

    @constraint(model, P * x[1:nc] .>= f_u)
    @constraint(model, P * x[nc+1:end] .<= f_l)

    @constraint(model, x[1:nc] .>= x[nc+1:end])

    for (i, ci) in enumerate(eachcol(C)), (j, cj) in enumerate(eachcol(C))
        if sum(ci .< cj) == 1 && sum(ci .== cj) == (size(C, 1) - 1)
            @constraint(model, x[i] <= x[j])
            @constraint(model, x[i+nc] <= x[j+nc])
        end
    end

    JuMP.optimize!(model)

    w_u = value.(x[1:nc])
    w_l = value.(x[nc+1:end])
    return w_u, w_l
end
