module IntervalPredictorModel
# ==============================================================================
using Convex
using SCS
# ==============================================================================
using ..SurvivalSignatureUtils
using ..Structures: Points, Method, PredictorModel
using ..BasisFunction

# needed for monotonicity_constraints
include("../src/rbf/radialbasisfunctions.jl")

# ==============================================================================
export intervalPredictor
# ==============================================================================

# no idea what is happening here at the momment
function intervalPredictor(
    evaluated_points::Points,
    upper_bound::AbstractArray,
    lower_bound::AbstractArray,
    centers::AbstractArray,
    shape_parameter::Union{Number,AbstractArray},
    weights::AbstractArray,
    method::Method,
)
    basis = BasisFunction.basis(
        method.basis_function_method,
        shape_parameter,
        evaluated_points.coordinates,
        centers,
        method.smoothness_factor,
    )

    num_centers = size(centers, 2)
    x = Variable(num_centers)               # requires Convex
    y = Variable(num_centers)               # requires Convex

    con = monotonicity_constraints(centers)

    con_x = x[con[1, :]] <= x[con[2, :]]
    con_y = y[con[1, :]] <= y[con[2, :]]

    con_ipm = x >= y
    con_u = (basis * x) >= upper_bound
    con_l = (basis * y) <= lower_bound

    N = size(basis, 2)

    problem = minimize(sum(basis * (x - y)) / N, [con_x, con_y, con_ipm, con_u, con_l])

    Convex.solve!(problem, SCS.Optimizer; silent=true)

    w_u = Convex.evaluate(x)
    w_l = Convex.evaluate(y)

    return PredictorModel(evaluated_points, centers, shape_parameter, weights, w_u, w_l)   #struct
end

# ==============================================================================
end