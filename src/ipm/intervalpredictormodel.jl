struct IntervalPredictorModel
    c::AbstractMatrix
    σ::Vector{<:Real}
    wmax::Vector{<:Real}
    wmin::Vector{<:Real}
end

function IntervalPredictorModel(
    X::AbstractMatrix,
    f_u::AbstractVector,
    f_l::AbstractVector,
    centers::AbstractMatrix,
    σ::AbstractVector,
)
    P = basis(X, centers, σ)

    nc = size(centers, 2)
    x = Variable(nc)
    y = Variable(nc)

    con = monotonicity_constraints(centers)

    con_x = x[con[1, :]] <= x[con[2, :]]
    con_y = y[con[1, :]] <= y[con[2, :]]

    con_ipm = x >= y
    con_u = (P * x) >= f_u
    con_l = (P * y) <= f_l

    n = size(P, 2)

    problem = minimize(sum(P * (x - y)) / n, [con_x, con_y, con_ipm, con_u, con_l])

    solve!(problem, SCS.Optimizer; silent_solver=true)

    w_u = Convex.evaluate(x)
    w_l = Convex.evaluate(y)

    return IntervalPredictorModel(centers, σ, w_u, w_l)
end

function evaluate(ipm::IntervalPredictorModel, x::AbstractVector)
    Ψ = basis(x, ipm.c, ipm.σ)

    ub = (Ψ * ipm.wmax)[1]
    lb = (Ψ * ipm.wmin)[1]

    return lb, ub
end

function evaluate(ipm::IntervalPredictorModel, x::Tuple)
    return evaluate(ipm, [x...])
end
