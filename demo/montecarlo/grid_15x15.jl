using Distributed
using ClusterManagers
using JLD2 # extra dependency

addprocs(SlurmManager(48); N=2, exeflags="--project")

@everywhere begin
    using SurvivalSignature

    n = 15
    m = 15
    adj = gridnetwork(n, m)
    types = Dict(1 => collect(3:2:((n * m) - 1)), 2 => collect(2:2:((n * m) - 1)))
    φ = s_t_connectivity(
        [2:((n * m) - 1);], findall(!iszero, adj[1, :]), findall(!iszero, adj[:, n * m])
    )
end

N = 10000
covtol = 1e-3

Φ, cov = survivalsignature(adj, types, φ, N, covtol, percolation_preprocessor!)

jldsave("grid-network-15x15-MC-10000"; Φ, cov)
