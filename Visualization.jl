module Visualization
# the purpose of this module is to house the functions related to plotting and
# and text outputs
# ==============================================================================

using ..Structures: System, Simulation, Method

# ==============================================================================

export printDetails

# ============================= TEXT ===========================================

function removeDashAndCapitalize(str::String)
    return join([titlecase(word) for word in split(str, "-")], " ")
end

function printDetails(sys::System, sim::Simulation, method::Method)
    println("==================================================")
    println("   Mode: $(removeDashAndCapitalize(method.simulation_method,))")
    println("--------------------------------------------------")
    println("   Parameters")
    println("..................................................")
    println("\tSamples: $(sim.samples)")
    println("\tVariation Tolerance: $(sim.variation_tolerance)")
    if method.simulation_method == "monte-carlo"
        nothing
    else
        println("\tWeight Change Tolerance: $(sim.weight_change_tolerance)")
        if method.simulation_method == "interval-predictor"
            println("\tConfidence Interval: $(sim.confidence_interval)")
        end
    end
    println("--------------------------------------------------")

    if method.simulation_method == "monte-carlo"
        nothing # monte-carlo doesnt use these methods
    else
        println("   Methods:")
        println("..................................................")
        println(
            "\tShape Parameter: $(removeDashAndCapitalize(method.shape_parameter_method))"
        )

        println(
            "\tBasis Function: $(removeDashAndCapitalize(method.basis_function_method))"
        )

        if method.basis_function_method == "matern"
            println("\t\tSmoothness Factor: $(method.smoothness_factor)")
        end

        println(
            "\tStarting Points: $(removeDashAndCapitalize(method.starting_points_method))"
        )
        println("\tCenter Points: $(removeDashAndCapitalize(method.centers_method))")

        println("\tWeight Change: $(removeDashAndCapitalize(method.weight_change_method))")
        println("==================================================")
        println("")
    end

    return nothing
end

# ============================ PLOTS ===========================================

function plotError()
    # plot errors of multiple systems (with various shape parameters)
    return nothing
end

# ==============================================================================

end