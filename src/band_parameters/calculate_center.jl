# band_parameters/center.jl

"
    calculate_center(
        contrem_spectrum::AbstractVector{T},
        λ::AbstractVector{T},
        λmin::T,
        λmax::T
    )::T where T <: AbstractFloat

Calculates the center of a spectral absorption band.

# Arguments
- `spectrum::AbstractVector{T}`: Spectral data containing an absorption.
- `λ::AbstractVector{T}`: All wavelengths.
- `λmin::T`: Minimum wavelength (matching units of spectrum) to use in
             calculation.
- `λmax::T`: Maximum wavelength to use in calculation.

# Returns
- `center::Float64`: Absorption Center.
- `xfit::Vector{Float64}`: X-Data used for fit.
- `yfit::Vector{Float64}`: Fitted Y-Data.
"
function calculate_center(
    contrem_spectrum::AbstractVector{T},
    λ::AbstractVector{T},
    λmin::V,
    λmax::V
)::Tuple{Float64, Vector{Float64}, Vector{Float64}} where {
    T<:AbstractFloat,
    V<:Real
}
    # Finding wavelength indices and real wavelengths from guesses.
    λmin_index, λmin = findλ(λ, λmin)
    λmax_index, λmax = findλ(λ, λmax)
    λindices = λmin_index:1:λmax_index

    N = 4  # Fit order

    absorption_spec = contrem_spectrum[λindices]
    absorption_λ = λ[λindices]
    β = polyfit(absorption_λ, absorption_spec, N)
    xfit = LinRange(λmin, λmax, 100) |> collect
    X = vandermonde(xfit, N)
    yfit = X * β
    center = xfit[argmin(yfit)]
    if center <= λmin
        return NaN64, xfit, yfit
    else
        return center, xfit, yfit
    end
end


function calculate_center(
    contrem_spec_cube::AbstractArray{T},
    λ::AbstractVector{T},
    λmin::V,
    λmax::V
)::Tuple{Matrix{Float64}, Vector{Float64}, Array{Float64}} where {
    T<:AbstractFloat,
    V<:Real
}
    
    # Finding wavelength indices and real wavelengths from guesses.
    λmin_index, λmin = findλ(λ, λmin)
    λmax_index, λmax = findλ(λ, λmax)
    λindices = λmin_index:1:λmax_index

    N = 4  # Fit order

    n, m, _ = size(contrem_spec_cube)

    absorption_spec_cube = contrem_spec_cube[:, :, λindices]
    absorption_λ = λ[λindices]

    β = polyfit(absorption_spec_cube, absorption_λ, N)

    xfit = LinRange(λmin, λmax, 100) |> collect
    X = vandermonde(xfit, N)
    yfit = X * β'

    yfit_cube = permutedims(reshape(yfit, size(xfit, 1), n, m), (2, 3, 1))
    center_idx = dropdims(mapslices(argmin, yfit_cube; dims=3); dims=3)
    center_map = xfit[center_idx]
    center_map[center_map.==λmin] .= NaN64
    center_map[center_map.==λmax] .= NaN64

    return center_map, xfit, yfit_cube
end


function calculate_center(
    spec::Spectrum{T},
    λmin::V,
    λmax::V
)::Tuple{T, AbstractVector{T}, AbstractVector{T}} where {
    T<:AbstractFloat,
    V<:Real
}
    return calculate_center(
        spec.contrem,
        spec.λ,
        λmin,
        λmax
    )
end

function calculate_center(
    spec_cube::SpectralCube{T},
    λmin::V,
    λmax::V
)::Tuple{Matrix{Float64}, Vector{Float64}, Array{Float64, 3}} where {
    T<:AbstractFloat,
    V<:Real
}
    return calculate_center(
        spec_cube.contrem,
        spec_cube.λ,
        λmin,
        λmax
    )
end


