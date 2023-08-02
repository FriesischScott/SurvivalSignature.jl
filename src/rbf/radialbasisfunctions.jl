function basis(
    X::Union{AbstractMatrix,AbstractVector}, Y::AbstractMatrix, σ::AbstractVector
)
    Ψ =
        exp.(
            -[
                sum((x .- c) .^ 2 ./ 2σ .^ 2) for
                (x, c) in Iterators.product(eachcol(X), eachcol(Y))
            ]
        )
    return Ψ ./ sum(Ψ; dims=2)
end

# function lsqr(P::AbstractMatrix, f::AbstractVector, con::Matrix{Int})
#     model = Model(SCS.Optimizer)
#     set_silent(model)
#     @variable(model, x[1:size(P, 2)])
#     @variable(model, res[1:length(f)])
#     @constraint(model, res .== P * x - f)

#     @constraint(model, x[con[1, :]] .<= x[con[2, :]])

#     @objective(model, Min, sum(res .^ 2))

#     JuMP.optimize!(model)
#     return value.(x)
# end

function lsqr(P::AbstractMatrix, f::AbstractVector, con::Matrix{Int})
    x = Variable(size(P, 2))

    con = (x[con[1, :]] - x[con[2, :]]) <= 0

    problem = minimize(sumsquares(P * x - f), con)
    solve!(problem, SCS.Optimizer; silent_solver=true)

    return Convex.evaluate(x)
end

function monotonicity_constraints(centers::AbstractMatrix)
    con = map(
        Iterators.product(enumerate(eachcol(centers)), 1:size(centers, 1))
    ) do ((i, ci), k)
        d = Not(k)
        f = findfirst(c -> ci[k] < c[k] && ci[d] == c[d], eachcol(centers))
        return [i, f]
    end

    con = filter(c -> !any(isnothing.(c)), con)
    return hcat(con...)
end
