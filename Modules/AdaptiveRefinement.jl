module AdaptiveRefinement
# ==============================================================================

using ForwardDiff
using NearestNeighbors
using InvertedIndices
using ProgressMeter

# ==============================================================================

using ..SurvivalSignatureUtils
using ..Structures: Points, System, Simulation, Methods
using ..Evaluation
using ..BasisFunction
using ..ShapeParameter
using ..Error

# access to monotonicity_constraints and lsqr 
include("../src/rbf/radialbasisfunctions.jl")

# ==============================================================================

export adaptiveRefinement

# ==============================================================================

function adaptiveRefinement(
    total_points::Points,
    evaluated_points::Points,
    sys::System,
    sim::Simulation,
    method::Methods,
    weights::AbstractArray,
    centers::AbstractArray,
    constraints::AbstractArray,
    shape_parameter::Union{Number,AbstractArray},
)
    total_basis, _ = BasisFunction.basis(
        method.basis_function_method, shape_parameter, total_points.coordinates, centers
    )

    x = evaluated_points # evaluated points
    l_max = maximumLength(total_points)
    candidates, cand_idx = remainingCandidates(total_points, x)

    prog = ProgressMeter.ProgressThresh(
        sim.weight_change_tolerance, "\tAdaptive Refinement:"
    )
    stop = 0
    while stop < 2
        function s(a::Points)
            return (BasisFunction.basis(method.basis_function_method, shape_parameter, a.coordinates, centers)[1] * weights)[1]
        end

        # exploration score
        idx, D = explorationScore(x, candidates)

        # exploitation score
        R = exploitationScore(total_basis, x, s, weights, candidates, cand_idx, idx)

        # weight function
        w = weightFunction(idx, l_max)

        # hybrid score function
        J = hybridScore(D, R, w)

        optimal_candidate, optimal_idx = optimalPoint(total_points, J, candidates)

        true_value, coefficient_variation = Evaluation.computeSurvivalSignatureEntry(
            sys, sim, optimal_candidate
        )

        # update the computed state_vectors
        x.coordinates = hcat(x.coordinates, optimal_candidate)

        x.solution = vcat(x.solution, true_value)
        x.confidence = vcat(x.confidence, coefficient_variation)
        x.idx = vcat(x.idx, optimal_idx)

        candidates, cand_idx = remainingCandidates(total_points, x)

        old_weights = weights

        # recompute the basis Function
        basis, _ = BasisFunction.basis(
            method.basis_function_method, shape_parameter, x.coordinates, centers
        )

        # update weights
        weights = lsqr(basis, x.solution, constraints)

        if Error.calculateError(method.weight_change_method, weights, old_weights) <
            sim.weight_change_tolerance
            stop += 1
        else
            stop = 0
        end

        stop != 1 && ProgressMeter.update!(
            prog,
            Error.calculateError(method.weight_change_method, weights, old_weights),
        )
    end

    upper_bound = min.(x.solution .+ (x.solution .* x.confidence), 1.0)
    lower_bound = max.(x.solution .- (x.solution .* x.confidence), 0.0)

    replace!(upper_bound, NaN => 1 / (sim.samples + 1))
    replace!(lower_bound, NaN => 0.0)

    return x, weights, upper_bound, lower_bound
end

function explorationScore(x::Points, candidates::AbstractArray)
    tree = NearestNeighbors.KDTree(x.coordinates)
    idx, dist = NearestNeighbors.nn(tree, candidates)
    return idx, dist
end

function exploitationScore(
    total_basis::AbstractArray,
    X::Points,
    func::Function,
    weights::AbstractArray,
    candidates::AbstractArray,
    cand_idx::InvertedIndex,
    nearest_neighbor_idx::AbstractArray,
)
    # Combine coordinates and solutions for local function s(a)
    combined = [(X.coordinates[:, i], X.solution[i]) for i in 1:size(X.coordinates, 2)]

    # Initialize gradient as a zeros array
    ∇s = [zeros(size(X.coordinates, 1)) for _ in 1:size(X.coordinates, 2)]
    ∇s = map(ForwardDiff.DiffResults.GradientResult, ∇s)  # Convert to gradient object

    # Compute gradients for each set of coordinates in `combined` with respect to `func`.
    # For each pair of coordinates (`coords`) and corresponding solution (`sol`),
    # a new `Points` instance is created with `coords` as `y` and `sol` as `solution`.
    # The gradient is calculated and stored in the preallocated gradient result `r`.
    ∇s = map(
        ((coords, sol), r) ->
            ForwardDiff.gradient!(r, y -> func(Points(y, nothing, sol, nothing)), coords),
        combined,
        ∇s,
    )

    gradient_values = ForwardDiff.DiffResults.value.(∇s)
    ∇s = ForwardDiff.DiffResults.gradient.(∇s)

    # Exploration score function
    idx = nearest_neighbor_idx
    t = @views gradient_values[idx] .+ map(
        (∇, a) -> dot(∇, a), ∇s[idx], eachcol(candidates .- X.coordinates[:, idx])
    )

    return @views abs.(total_basis[cand_idx, :] * weights .- t)
end

function weightFunction(nearest_neighbor_distance::AbstractArray, l_max::Number)
    return (1 .- nearest_neighbor_distance ./ l_max)
end

function maximumLength(Ω::Points)
    lb = minimum(Ω.coordinates; dims=2)
    ub = maximum(Ω.coordinates; dims=2)

    return sqrt(sum((ub .- lb) .^ 2))
end

function hybridScore(
    exploration_score::AbstractArray,
    exploitation_score::AbstractArray,
    weight_function::AbstractArray,
)

    # normalization
    exploration_score = exploration_score ./ maximum(exploration_score)
    exploitation_score = exploitation_score ./ maximum(exploitation_score)

    return exploration_score + weight_function .* exploitation_score
end

function optimalPoint(Ω::Points, scores::AbstractArray, candidates::AbstractArray)
    _, idx = findmax(scores)
    optimal_candidate = candidates[:, idx]  # optimal candidate (state_vector)

    # convert from the candidates index to the equivalent Ω index
    optimal_candidate_idx = findfirst(
        x -> all(x .== optimal_candidate), eachcol(Ω.coordinates)
    )

    return optimal_candidate, optimal_candidate_idx
end

function remainingCandidates(Ω::Points, x::Points)

    # remaining values
    cand_idx = InvertedIndices.Not(x.idx)
    candidates = Ω.coordinates[:, cand_idx]

    return candidates, cand_idx
end

function remainingCandidates(candidates::AbstractArray, remove_points::AbstractArray)
    candidates = candidates[.!in.(1:size(candidates, 1), Ref(remove_points)), :]
    return candidates, cand_idx
end

# ==============================================================================

end
