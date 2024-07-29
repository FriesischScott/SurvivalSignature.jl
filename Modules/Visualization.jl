module Visualization
# the purpose of this module is to house the functions related to plotting and
# and text outputs

# ==============================================================================

using Plots

# ==============================================================================

using ..Structures: System, Simulation, Methods

using ..Structures: SimulationType, MonteCarloSimulation, IntervalPredictorSimulation
using ..Structures: CentersMethod, GridCenters

# ==============================================================================

export printDetails, plotError

# ============================= TEXT ===========================================

# needs to be completly rewritten for refactorization
# make multiple for each simulation type
function printDetails(sys::System, sim::Simulation, method::Methods)
    println("==================================================")
    println("   Mode: $(nameof(typeof(method.simulation_method)))")
    println("--------------------------------------------------")
    println("   Parameters")
    println("..................................................")
    println("\tSamples: $(sim.samples)")
    println("\tVariation Tolerance: $(sim.variation_tolerance)")

    println("--------------------------------------------------")

    if typeof(method.simulation_method) == MonteCarloSimulation()
        nothing # monte-carlo doesnt use these methods
    else
        println("   Methods:")
        println("..................................................")
        println("\tShape Parameter: $(nameof(typeof(method.shape_parameter_method)))")

        println("\tBasis Function: $(nameof(typeof(method.basis_function_method)))")

        println("\tStarting Points: $(nameof(typeof(method.starting_points_method)))")
        println("\tCenter Points: $(nameof(typeof(method.centers_method)))")
        if typeof(method.centers_method) == GridCenters()
            println("\t\tCenters Interval: $(method.centers_method.confidence_interval)")
        end

        println("\tWeight Change: $(nameof(typeof(method.weight_change_method)))")
        println("==================================================")
        println("")
    end

    return nothing
end

# ============================ PLOTS ===========================================

function plotErrorComparison(errors::Vararg{Dict{String,Union{Float64,String}}})
    names = [get(e, "name", "Unknown") for e in errors]
    x = unique([key for e in errors for key in keys(e) if key != "name"])

    plt = plot()  # Initialize a new plot

    for error_type in x
        plot_data = []
        plot_names = []

        for (e, name) in zip(errors, names)
            if haskey(e, error_type)
                push!(plot_data, e[error_type])
                push!(plot_names, name)
            end
        end

        scatter!(
            plot_names,
            plot_data;
            label=error_type,
            ylabel="Error",
            legend=true,
            xticks=:auto,
            xrotation=45,
        )
    end

    return plt
end

function plotErrorSequence(errors::Vararg{Vector{Dict{String,Union{Float64,String}}}})
    names = [get(e, "name", "Unknown") for e in errors]
    x = unique([key for e in errors for key in keys(e) if key != "name"])

    plt = plot()  # Initialize a new plot

    for error_type in x
        plot_data = []
        plot_names = []

        for (e, name) in zip(errors, names)
            if haskey(e, error_type)
                push!(plot_data, e[error_type])
                push!(plot_names, name)
            end
        end

        scatter!(
            plot_names,
            plot_data;
            label=error_type,
            ylabel="Error",
            legend=true,
            xticks=:auto,
            xrotation=45,
        )
    end

    return plt
end

function plotError(errors::Dict{String,Float64}...)

    # Extract unique error types
    error_types = unique([keys(d) for d in errors]...)

    # Determine number of rows and columns for subplots
    num_plots = length(error_types)
    num_cols = ceil(Int, sqrt(num_plots))
    num_rows = ceil(Int, num_plots / num_cols)

    # Create a plot array to hold individual plots
    plot_array = []

    # Loop over each error type
    for err_type in error_types
        # Extract data for the current error type from each dictionary
        plot_data = [get(d, err_type, 0) for d in errors]

        # Create a new plot for the current error type
        p = bar(plot_data; label=err_type)
        push!(plot_array, p)
    end

    # Plot all subplots in a grid layout
    return plot(plot_array...; layout=(num_rows, num_cols), legend=:topleft)
end

function plotWeightChange(weight_change::Vector...)
    return nothing
end
# ==============================================================================

end