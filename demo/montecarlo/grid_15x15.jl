using Distributed
using ClusterManagers
using JLD2 # extra dependency

addprocs(SlurmManager(48); N=2, exeflags="--project")

@everywhere begin
    using SurvivalSignature

    n = 15
    m = 15
    adj = gridnetwork(n, m)
    types = Dict(1 => collect(1:2:(n * m)), 2 => collect(2:2:(n * m)))
    φ = s_t_connectivity([1:(n * m);], [1], [n * m])
end

N = 10000
covtol = 1e-3

Φ, cov = survivalsignature(adj, types, φ, N, covtol, percolation_preprocessor!)

jldsave("grid-network-15x15-MC-10000.jld2"; Φ, cov)
