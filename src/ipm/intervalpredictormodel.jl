struct IntervalPredictorModel
    c::AbstractMatrix
    σ::Vector{<:Real}
    wmax::Vector{<:Real}
    wmin::Vector{<:Real}
end

function IntervalPredictorModel(X::AbstractMatrix, f_u::AbstractVector, f_l::AbstractVector, centers::AbstractMatrix, σ::AbstractVector)
    P = basis(X, centers, σ)

    nc = size(centers, 2)

    model = Model(SCS.Optimizer)
    set_silent(model)
    @variable(model, x[1:nc*2])
    @variable(model, spread[1:length(f_u)])

    @constraint(model, spread .== P * (x[1:nc] - x[nc+1:end]))
    @objective(model, Min, mean(spread))

    @constraint(model, P * x[1:nc] .>= f_u)
    @constraint(model, P * x[nc+1:end] .<= f_l)

    @constraint(model, x[1:nc] .>= x[nc+1:end])

    λ = size(centers, 1) - 1

    for (i, ci) in enumerate(eachcol(centers)), (j, cj) in enumerate(eachcol(centers))
        if sum(ci .< cj) == 1 && sum(ci .== cj) == λ
            @constraint(model, x[i] <= x[j])
            @constraint(model, x[i+nc] <= x[j+nc])
        end
    end

    JuMP.optimize!(model)

    w_u = value.(x[1:nc])
    w_l = value.(x[nc+1:end])

    return IntervalPredictorModel(centers, σ, w_u, w_l)
end

function evaluate(ipm::IntervalPredictorModel, x::AbstractVector)
    Ψ = basis(x, ipm.c, ipm.σ)

    ub = (Ψ*ipm.wmax)[1]
    lb = (Ψ*ipm.wmin)[1]

    return lb, ub
end

function evaluate(ipm::IntervalPredictorModel, x::Tuple)
    return evaluate(ipm, [x...])
end
