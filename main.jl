using SurvivalSignature
using ProgressMeter
using JLD2              # needed for load
# ==============================================================================

include("Modules/Import.jl")
using .Import

using ..SurvivalSignatureUtils
using ..Structures: System, Simulation, Model, Methods, StringMethods
using ..StructureCompilation
using ..Error: calculateError
using ..Systems
using ..Visualization: plotError

function main()
    # ================================ INPUTS ==================================

    system_type::String = "grid"                 # ["grid"]
    percolation::Bool = true

    # currently issues when n and m are greater than 11.  - using behrensdorf
    # also an error if n and m are too small (less than 7 each)
    n::Int = 15
    m::Int = 15

    # Simulation Parameters
    samples::Union{Int,Vector{Int}} = [1, 10, 100, 1000]
    covtol::Float64 = 1e-3                       # coeficient of varriation tolerance
    wtol::Float64 = 1e-3                         # weight change tolerance
    ci::Vector{Int} = [15, 15]          # confidence interval 

    # Methods
    # ["monte-carlo", "radial-basis-function", "interval-predictor"]
    simulation_method::Union{String,Vector{String}} = "interval-predictor"
    starting_points_method::Union{String,Vector{String}} = "grid-aligned"   # ["grid-aligned"]
    centers_method::Union{String,Vector{String}} = "grid-aligned"           # ["grid-aligned"]
    weight_change_method::Union{String,Vector{String}} = "norm"             # ["norm"]
    # only 'hardy' and 'rippa' seems to function correctly (i.e. result in reasonable error)
    #["hardy", "franke", "kuo", "rippa", "behrensdorf"]
    shape_parameter_method::Union{String,Vector{String}} = ["hardy", "behrensdorf"]
    basis_function_method::Union{String,Vector{String}} = "gaussian"       #["gaussian", "matern", "behrensdorf", "mq"]
    smoothness_factor::Int = 3                                             # only requred for matern [1, 2, 3]
    #                                                                     # and mq [1, 2]

    # Error
    error_type::Union{Vector{String},String} = ["rmse", "rae"]

    # ======================== STRUCT REFINEMENT ===============================

    sys::System = Systems.generateSystem(system_type, n, m; percolation_bool=percolation)

    sims::Union{Vector{Simulation},Simulation} = StructureCompilation.compileSimulation(
        samples, covtol, wtol, ci
    )

    str_methods::Union{Vector{StringMethods},StringMethods} = StructureCompilation.compileMethods(
        simulation_method,
        starting_points_method,
        centers_method,
        weight_change_method,
        shape_parameter_method,
        basis_function_method,
        smoothness_factor,
    )

    # =============================== SIMULATE =================================

    signatures::Union{Vector{Model},Model} = Simulate.simulate(
        simulation_method, sys, sims, str_methods
    )

    # ============================= "TRUE SOLUTION" ============================

    # sim_mc = Structures.Simulation(samples, covtol, wtol, ci, nothing) 
    # str_method_mc = Structures.StringMethods("monte-carlo")
    # signature_mc = Simulate.simulate("monte-carlo", sys, sim_mc, str_method_mc)

    #true values - apparently this outputs a Φ - which is the true value solutions
    @load "demo/data/grid-network-15x15-MC-10000.jld2"

    # =============================== ERROR ====================================

    println("Calculating Error...")
    println("--------------------------------------------------------")

    signatures, errors = Error.calculateError(error_type, signatures, Φ)

    println("--------------------------------------------------------")

    #plt = Visualization.plotErrorComparison(errors...)

    #display(plt)

    println("Errors Calculated.")
    return println("")
end

# =============================== RUN ==========================================
main()

#
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