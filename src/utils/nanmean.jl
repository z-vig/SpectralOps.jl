# utils/nanmath.jl

"
    nanmean(A, dim)

Takes the mean of A along dim, excluding NaN values.
"
function nanmean(A, dim)
    isnanmask = .!isnan.(A)
    sums = sum(x -> isnan(x) ? 0.0 : x, A; dims=dim)
    counts = sum(isnanmask; dims=dim)
    return sums ./ counts
end