function distance(X::AbstractMatrix, Y::AbstractMatrix, σ::Vector)
    return [sqrt(sum((x .- c) .^ 2 ./ 2σ .^ 2)) for (x, c) in Iterators.product(eachcol(X), eachcol(Y))]
end

function gaussian(x::Real)
    return [exp(-x^2)]
end
