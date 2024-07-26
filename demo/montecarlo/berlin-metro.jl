using Distributed
using Distributions
using JLD2

addprocs(12; exeflags="--project")

@everywhere begin
    using SurvivalSignature

    include("../models/berlin-metro/berlin-metro-three-types.jl")

    φ = efficiency(nodes, A, 0.5)
end

N = 10000
covtol = 1e-3

Φ, cov = survivalsignature(A, types, φ, N, covtol, percolation_preprocessor!)

jldsave("berlin-metro-three-types-0-5-10000-MC.jld2"; Φ, cov)
