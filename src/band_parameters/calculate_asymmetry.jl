# band_parameters/asymmetry.jl

"
    calculate_asymmetry(
        spectrum::AbstractVector{T},
        λmin::T,
        λmax::T
    )::T where T <: AbstractFloat

Calculates the asymmetry of a spectral absorption band.

# Arguments
- `spectrum::AbstractVector{T}`: Spectral data containing an absorption.
- `λmin::T`: Minimum wavelength (matching units of spectrum) to use in
             calculation.
- `λmax::T`: Maximum wavelength to use in calculation.
"
function calculate_asymmetry(
    spectrum::AbstractVector{T},
    λmin::V,
    λmax::V
)::T where {
    T <: AbstractFloat,
    V <: Real
}
    @warn "The `calculate_asymmetry` function has not been implemented yet, "*
          "returning NaN64"
    return NaN
end
