# smoothing/outlier_removal.jl

"""
    remove_outliers(
        spectrum::AbstractVector{T};
        threshold=2
    ) :: Vector{Float64} where T<:AbstractFloat

Removes outliers from a spectrum using the threshold method. The function
returns the spectrum with the outliers removed.

# Arguments
- `spectrum::Vector{<:AbstractFloat}`: The input spectrum with outliers.
- `threshold::Float16`: Number of standard deviations to use as outlier
  threshold.
"""
function remove_outliers(
    spectrum::AbstractVector{T};
    threshold::V = 2
)::Vector{Float64} where {
    T<:AbstractFloat,
    V<:Real
}
    bs = max(round_to_odd(length(spectrum)*0.25), 3)
    if iseven(bs) bs -= 1 end
    μ, σ, idx = moving_average(
        spectrum,
        box_size=bs,
        edge_handling="extrapolate",
        rm_outliers=false
    )

    zscore = (spectrum .- μ) ./ σ
    outlier_idx = abs.(zscore) .> threshold

    neighbors = stack([
        circshift!(spectrum, -1),
        circshift!(spectrum, 1)
    ], dims=2)

    neighbors[1, 2] = NaN
    neighbors[end, 1] = NaN

    replacement = nanmean(neighbors, 2)
    
    spectrum[outlier_idx] .= replacement[outlier_idx]
    return spectrum
end

function remove_outliers(
    cube::Array{T, 3},
    threshold::V = 2
)::Array{Float64, 3} where {
    T<:AbstractFloat,
    V<:Real
}
    ax = axes(cube)
    no_outliers = map(CartesianIndices(ax[1:2])) do I
        x,y = Tuple(I)
        spec = cube[x,y,:]
        return remove_outliers(spec, threshold=threshold)
    end

    return make3d(no_outliers)
end