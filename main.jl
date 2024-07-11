using SurvivalSignature
using ProgressMeter
using JLD2              # needed for load
using Statistics        # needed for mean

using Profile
using ProfileView
# ==============================================================================

include("Modules/Import.jl")
using .Import

using ..SurvivalSignatureUtils
using ..Structures: System, Simulation, Method, Model
using ..BasisFunction
using ..Error: calculateError
using ..Systems
# ================================ INPUTS ======================================

system_type = "grid"                 # ["grid"]
percolation_bool = true

# currently issues when n and m are greater than 11.  - using behernsdorf
# also an error if n and m are too small (less than 7 each)
n = 15
m = 15

# Simulation Parameters
samples = 10^2
covtol = 1e-3   # coeficient of varriation tolerance
wtol = 1e-3    # weight change tolerance
ci = [15, 15]   # confidence interval 

# Methods
simulation_method = "interval-predictor"  # ["monte-carlo", "radial-basis-function", "interval-predictor"]
starting_points_method = "grid-aligned"   # ["grid-aligned"]
centers_method = "grid-aligned"           # ["grid-aligned"]
weight_change_method = "norm"             # ["norm"]
shape_parameter_method = "behrensdorf"        #["hardy", "franke", "kuo", "rippa", "behrensdorf"]
basis_function_method = "matern"          #["gaussian", "matern", "behrensdorf"]
smoothness_factor = 3                    # only requred for matern [1, 2, 3]

# ======================== STRUCT REFINEMENT ===================================

# 'behrensdorf' has a different basis function method, and must be adjusted accordingly
if lowercase(shape_parameter_method) == "behrensdorf"
    basis_function_method = "behrensdorf"
end

sys = Systems.generateSystem(system_type, n, m; percolation_bool=percolation_bool)
sim = Structures.Simulation(samples, covtol, wtol, ci, nothing)
method = Structures.Method(
    simulation_method,
    starting_points_method,
    centers_method,
    weight_change_method,
    shape_parameter_method,
    basis_function_method,
    smoothness_factor,
)

# ============================= SIMULATE =======================================

# this is slower than behernsdorfs version 
# possibly just due to more if statements because of methods, 
# but i cant get the Profiler to work
signature = Simulate.simulate(sys, sim, method)

method = Structures.Method("monte-carlo")
# signature_mc = Simulate.simulate(sys, sim, method)

# =============================== ERROR ========================================

# true values - apparently this outputs a Φ - which is the true value solutions
@load "demo/data/grid-network-15x15-MC-10000.jld2"

function behrensdorfRSME(signature::Model, comparison::Matrix) #rmse()
    w = signature.model.weights
    threshold = signature.sim.threshold

    e = Float64[]

    for idx in CartesianIndices(comparison)
        if sum([Tuple(idx)...] .- 1) < threshold
            continue
        end
        s = (BasisFunction.basis(signature.method.basis_function_method, signature.model.shape_parameter, [Tuple(idx)...], signature.model.centers, signature.method.smoothness_factor) * w)[1]

        push!(e, s - comparison[idx])
    end
    return sqrt(Statistics.mean(e .^ 2))
end
function behrensdorfRAE(signature::Model, comparison::Matrix)
    w = signature.model.weights
    threshold = signature.sim.threshold

    e = Float64[]
    y = Float64[]

    for idx in CartesianIndices(comparison)
        if sum([Tuple(idx)...] .- 1) < threshold
            continue
        end
        s = (BasisFunction.basis(signature.method.basis_function_method, signature.model.shape_parameter, [Tuple(idx)...], signature.model.centers, signature.method.smoothness_factor) * w)[1]

        push!(e, s - comparison[idx])
        push!(y, comparison[idx])
    end
    return Statistics.mean(abs.(e)) / Statistics.mean(abs.(y .- Statistics.mean(y)))
end

println("Calculating Error...")

error_val_rmse = Error.calculateError("rmse", signature.Phi.solution, Φ)
error_val_rae = Error.calculateError("rae", signature.Phi.solution, Φ)

error_behrensdorf_rmse = behrensdorfRSME(signature, Φ)
error_behrensdorf_rae = behrensdorfRAE(signature, Φ)

println("\trmse: $error_val_rmse")
println("\trae: $error_val_rae")
println("\tbehrensdorf_rmse: $error_behrensdorf_rmse")
println("\tbehrensdorf_rae: $error_behrensdorf_rae")

println("Errors Calculated.")
println("")

#
#
#
#
#
#
#
#
#
# TO DO:
#    Test current progress with all shape parameter methods
#                                        i.e. not just behrensdorf
#
#    make a simuate function which encompasses main.jl
#
#   make visualization functions 
#   make Statistics text print-out functions
#   make error modules
#   make monte carlo and radial basis functions
#   add more shape parameter calculations
#   
#
#   potentially make more weight change functions
#       possible overlap of the error module functions if designed correctly.
#
#  do other shape parameters need to be adjusted to work in interval predictors?

# convert state_vectors to include all state_vectors but with solution set to zero
#     remember remaining coordinates must be updated in accordance !!!
#
#
# FINISH: 
#
#   evaluateSurrogate function
#           use in AdaptiveRefinement function as well
#   best construction of Model struct
#   determine what w_u, w_l, f_u, and f_l are, and if they are used.
#     might be used for the equivalent to my evaluateSurrogate function
#     for the upper and lower bounds of the weights, but im not sure.
#     i cant seem to find those functions used.
#
#    speaking of: im not sure what this does:
#            w_u = Convex.evaluate(x)
#            w_l = Convex.evaluate(y)
#
#    x, y are just the 'Variable' of the number of centers. 
#
#    also they are the same thing, so shouldnt result in upper and lower bowers
#    despite how they are utilized.
#
#    f_u and f_l seem to be related to the solutions to get confidence intervals of some sort. 
#    possibly only needed for the interval predictors 
#    keep in mind when making the monteCarlo and radialBasisFunction simulations
#
#    use @examine to determine more about them.
#
#
#
# ISSUE WHEN THE GRID GETS TOO BIG IT SEEMS 
#