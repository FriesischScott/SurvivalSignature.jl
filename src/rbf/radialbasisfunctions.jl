function distance(X::AbstractMatrix, Y::AbstractMatrix, σ::AbstractVector)
    return [sqrt(sum((x .- c) .^ 2 ./ 2σ .^ 2)) for (x, c) in Iterators.product(eachcol(X), eachcol(Y))]
end

function distance(x::AbstractVector, Y::AbstractMatrix, σ::AbstractVector)
    return [sqrt(sum((x .- c) .^ 2 ./ 2σ .^ 2)) for c in eachcol(Y)]
end

function gaussian(x::Real)
    return exp(-x^2)
end

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
