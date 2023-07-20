module SurvivalSignature

using BigCombinatorics
using Distributed
using Distributions
using ForwardDiff
using InvertedIndices
using IterTools
using JuMP
using LinearAlgebra
using NearestNeighbors
using ProgressMeter
using Random
using SCS
using Statistics

export survivalsignature
export exactentry
export approximateentry
export gridnetwork
export percolation
export s_t_connectivity
export reliability
export percolation_preprocessor!

export IPMSurvivalSignature

include("ipm/intervalpredictormodel.jl")

struct IPMSurvivalSignature
    X::AbstractMatrix
    f::AbstractVector
    k::AbstractVector
    fc::Real
    ipm::IntervalPredictorModel
end

include("rbf/radialbasisfunctions.jl")

include("percolation.jl")
include("preprocessors.jl")
include("reliability.jl")
include("signature.jl")
include("structurefunctions.jl")
include("util.jl")

end # module
