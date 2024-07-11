module Import
# ==============================================================================

export SurvivalSignatureUtils,
    Structures,
    Visualization,
    Error,
    Systems,
    BasisFunction,
    ShapeParameter,
    StartingPoints,
    Centers,
    Evaluation,
    IntervalPredictorModel,
    AdaptiveRefinement,
    Simulate,
    monotonicity_constraints, #from final include
    lsqr                      #from final include
# ==============================================================================
include("SurvivalSignatureUtils.jl")
using .SurvivalSignatureUtils

include("Structures.jl")
using .Structures

include("Visualization.jl")
using .Visualization

include("BasisFunction.jl")
using .BasisFunction

include("Error.jl")
using .Error

include("Systems.jl")
using .Systems

include("ShapeParameter.jl")
using .ShapeParameter

include("StartingPoints.jl")
using .StartingPoints

include("Centers.jl")
using .Centers

include("Evaluation.jl")
using .Evaluation

include("IntervalPredictorModel.jl")
using .IntervalPredictorModel

include("AdaptiveRefinement.jl")
using .AdaptiveRefinement

include("Simulate.jl")
using .Simulate

# convert this to a Module: RadialBasisFunctions (?)
# needed for monotonicity_constraints and lsqr
include("../src/rbf/radialbasisfunctions.jl")
# ==============================================================================
end