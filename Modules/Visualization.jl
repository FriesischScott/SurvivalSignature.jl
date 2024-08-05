module Visualization
# the purpose of this module is to house the functions related to plotting and
# and text outputs

# ==============================================================================

using CairoMakie, LaTeXStrings

# ==============================================================================

using ..Structures: System, Simulation, Methods, Model

using ..Structures: SimulationType, MonteCarloSimulation, IntervalPredictorSimulation
using ..Structures: CentersMethod, GridCenters

using ..SurvivalSignatureUtils

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

    if typeof(method.simulation_method) == MonteCarloSimulation
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
function looped_variable(structs...)
    # Get the field names of the first struct
    fields = fieldnames(typeof(structs[1]))

    # Ensure all structs have the same fields
    for s in structs[2:end]
        if fieldnames(typeof(s)) != fields
            error("Structs have different fields")
        end
    end

    # Find the differing field
    differing_field = nothing
    for field in fields
        values = getfield.(structs, field)
        unique_values = unique(values)

        if length(unique_values) > 1
            if isnothing(differing_field)
                differing_field = (field, unique_values)
            else
                error("More than one differing field found")
            end
        end
    end

    if isnothing(differing_field)
        return false
    else
        return differing_field
    end
end

function extractError(
    errors::Union{Matrix{Dict{String,Float64}},Vector{Dict{String,Float64}}}
)
    # assumes all dicts have the same keys
    key_names = collect(Base.keys(errors[1]))

    extracted_errors = Vector{Tuple{String,Array{Float64}}}(undef, length(key_names))

    for (i, key) in enumerate(key_names)
        extracted_errors[i] = extractError(errors, key)
    end

    return extracted_errors
end

function extractError(errors::Vector{Dict{String,Float64}}, key::String)
    rows = length(errors)
    error_matrix = Array{Float64}(undef, rows)

    for i in 1:rows
        error_matrix[i] = get(errors[i], key, NaN)
    end

    return (key, error_matrix)
end

function extractError(errors::Matrix{Dict{String,Float64}}, key::String)
    rows, cols = size(errors)
    error_matrix = Array{Float64}(undef, rows, cols)

    for i in 1:rows
        for j in 1:cols
            error_matrix[i, j] = get(errors[i, j], key, NaN)
        end
    end

    return (key, error_matrix)
end

function removeDashAndCapitalize(str::String)
    return join([titlecase(word) for word in split(str, "-")], " ")
end

function removeUnderscoreAndCapitalize(str::String)
    return join([titlecase(word) for word in split(str, "_")], " ")
end

function latexSpaces(str::String)
    return replace(str, " " => "\\ ")
end

# ============================ PLOTS ===========================================

function plotError(
    axis::Makie.Axis,
    error::Vector{Float64},
    looped_variable::Vector{<:Number};
    labels::Union{String,Nothing}=nothing,
)
    return CairoMakie.lines!(
        axis,
        looped_variable,
        error;
        label=isnothing(labels) ? nothing : LaTeXStrings.latexstring(labels),
    )
end

function plotError(
    axis::Makie.Axis,
    error::Matrix{Float64},
    looped_variable::Vector{<:Number};
    labels::Union{Vector{String},Nothing}=nothing,
)
    for (i, row::Vector) in enumerate(eachrow(error))
        plotError(
            axis, row, looped_variable; labels=isnothing(labels) ? nothing : labels[i]
        )
    end
end

function plotError(
    signatures::Union{Model,Array{Model}}, errors::Array{Dict{String,Float64}}
)
    sim_variable = looped_variable(getfield.(signatures, :sim)...)
    method_variable = looped_variable(getfield.(signatures, :method)...)

    extracted_errors = extractError(errors)

    plt = Figure(; fonts=(; regular="CMU Serif"))

    for (i, (key, error)) in enumerate(extracted_errors)
        ax = Axis(
            plt[i, 1];
            title=LaTeXStrings.latexstring(key),
            xlabel=LaTeXStrings.latexstring(uppercase(string((sim_variable[1])))),
            xgridstyle=:dash,
            ygridstyle=:dash,
            xtickalign=1,
            xticksize=5,
            ytickalign=1,
            yticksize=5,
            xscale=log10,
            xminorticksvisible=true,
            yscale=log10,
            yminorticksvisible=true,
        )

        if typeof(method_variable) == Bool
            plotError(ax, error, sim_variable[2])
        else
            plotError(
                ax,
                error,
                sim_variable[2];
                labels=string.(nameof.(typeof.(method_variable[2]))),
            )
        end
        if typeof(method_variable) == Bool
        else
            Legend(plt[i, 2], ax; unique=true, merge=true)
        end
    end

    return plt
end
# ==============================================================================

end