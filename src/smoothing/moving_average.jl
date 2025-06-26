# smoothing/moving_average.jl

using DSP

"""
    moving_average(
        spectrum::Vector{<:AbstractFloat};
        box_size=5,
        edge_handling::String="extrapolate"
    )::Tuple{Vector{Float64},Vector{Float64},Vector{Int}}

Smooths the input spectrum using a moving average filter. The function returns
the smoothed spectrum, the standard deviation of the smoothed spectrum, and
the indices of the smoothed spectrum. The `box_size` parameter determines the
size of the moving average filter. The `edge_handling` parameter determines
how the edges of the spectrum are handled. The options are "mirror",
"extrapolate", "fill_ends", and "cut_ends". The "mirror" option mirrors the
spectrum at the edges, the "extrapolate" option extrapolates the spectrum at
the edges, the "fill_ends" option fills the ends of the spectrum with the
original data, and the "cut_ends" option cuts the ends of the spectrum.

# Arguments
- `spec::Vector{<:AbstractFloat}`: The input spectrum to be smoothed.
- `box_size::Int=5`: The size of the moving average filter.
- `edge_handling::String="extrapolate"`: The method used to handle the edges of
  the spectrum.
    - "mirror": Mirrors the spectrum at the edges.
    - "extrapolate": Extrapolates the spectrum at the edges.
    - "fill_ends": Fills the ends of the spectrum with the original data.
    - "cut_ends": Cuts the ends of the spectrum.

# Returns
- `μ`: Smoothed spectrum
- `σ`: Errors on the smoothed spectum
- `idx`: List of indices for the smoothed spectrum
"""
function moving_average(
    spectrum_original::Vector{<:AbstractFloat};
    box_size::Int=5,
    edge_handling::String="extrapolate",
    rm_outliers::Bool=true
)::Tuple{Vector{Float64},Vector{Float64},Vector{Int}}

    spectrum = copy(spectrum_original)

    if iseven(box_size)
        println("box_size must be even!")
    end

    box = ones(box_size)
    endcap = box_size ÷ 2

    if rm_outliers
        spectrum = remove_outliers(spectrum, threshold=1.5)
    end

    if edge_handling == "mirror"
        spectrum = [
            reverse(spectrum[begin:begin+(box_size-1)]);
            spectrum;
            reverse(spectrum[end-(box_size-1):end])
        ]
    elseif edge_handling == "extrapolate"
        # We are going to fix the number of points used for the linear
        # extrapolation based on the lenght of the spectrum (10% of the
        # spectrum length). Changing the box size will simply effect how many
        # points are extrapolated, not the number of points used for the
        # extrapolation.
        fit_order = 1
        edge_length = max(round(Int, length(spectrum)*0.1), 1)
        left_idx = 1:1+edge_length
        right_idx = length(spectrum)-edge_length:length(spectrum)
        
        fit_left = linearfit(
            Vector(left_idx),
            spectrum[left_idx],
            [i for i in -box_size+1:0]
        )
        fit_right = linearfit(
            Vector(right_idx),
            spectrum[right_idx],
            [i for i in length(spectrum):length(spectrum)+box_size-1]
        )

        spectrum = [
            fit_left;
            spectrum;
            fit_right
        ]
    end

    μ = conv(spectrum,box)[begin+endcap:end-endcap] ./ box_size
    μ² = conv(spectrum.^2,box)[begin+endcap:end-endcap] ./ box_size
    
    idx = 1:length(spectrum)

    #Removes convolution artifacts from the edges
    if edge_handling == "mirror" || edge_handling == "extrapolate"
        μ = μ[begin+box_size:end-box_size]
        μ² = μ²[begin+box_size:end-box_size]
        idx = idx[1:end-(2*box_size)]
    end

    σ = sqrt.(μ² .- μ.^2)

    if edge_handling == "fill_ends"
        μ[begin:begin+endcap] = spectrum[begin:begin+endcap]
        μ[end-endcap:end] = spectrum[end-endcap:end]
        σ[begin:begin+endcap] .= 0
        σ[end-endcap:end] .= 0
    elseif edge_handling == "cut_ends"
        μ = μ[begin+endcap:end-endcap]
        σ = σ[begin+endcap:end-endcap]
        idx = 1+endcap:length(spectrum)-endcap
    end

    return μ, σ, idx
end

function moving_average(
    cube::Array{T, 3},
    box_size::Int=5,
    edge_handling::String="extrapolate"
)::Tuple{Array{T, 3}, Array{T, 3}, Array{T, 3}} where T<:AbstractFloat

    ax = axes(cube)
    mvavg_res = map(CartesianIndices(ax[1:2])) do I
        x,y = Tuple(I)
        spec = cube[x,y,:]
        return moving_average(
            spec, box_size=box_size, edge_handling=edge_handling
        )
    end
    # Separate the matrix of tuples into three matrices
    μ, σ, idx = separate_tuple_matrix(mvavg_res)

    return μ, σ, idx
end