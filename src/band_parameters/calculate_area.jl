# band_parameters/area.jl

"
    calculate_area(
        spectrum::AbstractVector{T},
        λmin::T,
        λmax::T,
        spectral_resolution::T
    )::T where T <: AbstractFloat

Calculates the area of a spectral absorption band.

# Arguments
- `spectrum::AbstractVector{T}`: Spectral data containing an absorption.
- `λmin::T`: Minimum wavelength (matching units of spectrum) to use in
             calculation.
- `λmax::T`: Maximum wavelength to use in calculation.
- `spectral_resolution::T`: Spectral resolution in nm per band.

# Returns
- `area::Float64`: Absorption band area.
"
function calculate_area(
    contrem_spectrum::AbstractVector{T},
    λmin::V,
    λmax::V,
    spectral_resolution::T
)::Float64 where {
    T<:AbstractFloat,
    V<:Real
}
    # Finding wavelength indices and real wavelengths from guesses.
    λmin_index, λmin = findλ(λ, λrange[1])
    λmax_index, λmax = findλ(λ, λrange[2])
    λindices = λmin_index:1:λmax_index

    area_components = (1 .- contrem_spectrum[λindices]) .*
           spectral_resolution

    area = sum(area_components)

    return area
end

"
Convenience function to calculate band area using Spectrum object.
"
function calculate_area(
    spec::Spectrum{T},
    λmin::V,
    λmax::V
) where {
    T<:AbstractFloat,
    V<:Real
}
    return calculate_area(
        spec.contrem,
        λmin,
        λmax,
        spec.spec_res
    )
end
