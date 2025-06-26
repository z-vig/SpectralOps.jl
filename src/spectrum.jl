# Spectrum.jl

"
    Spectrum(
        data::V
        λ::V
        spectral_units::String
        λ_units::String
    ) {T<:AbstractFloat, V<:AbstractVector}

Represents a spectrum over a set of wavelengths.

# Fields
- `data::V`: Spectral data.
- `λ::V`: Wavelength data.
- `spectral_units::String`: Units of spectral data (One of: Reflectance, RFL,
                            Radiance or RDN).
- `λ_units::String`: Units of wavelength (one of: nanometers, nm, microns, um).
- `spec_res::T`: Spectral resolution of the spectrum.
- `no_outliers::V`: Spectrum but with outliers removed.
- `smoothed::V`: Spectrum with moving average applied.
- `smoothed_err::V`: 1σ Error on the moving average.
- `contrem::V`: Spectrum with the lunar double line continuum removed.
- `continuum::V`: Spectral continuum values.
"
struct Spectrum{T<:AbstractFloat, V<:AbstractVector{T}}
    data::V
    λ::V
    spectral_units::String
    λ_units::String
    spec_res:: T
    no_outliers::V
    smoothed::V
    smoothed_err::V
    contrem::V
    continuum::V

    function Spectrum(
        data::V,
        λ::V,
        spectral_units::String,
        λ_units::String
    ) where {
        T<:AbstractFloat,
        V<:AbstractVector{T}
    }
        validate_spectrum(data, λ, spectral_units, λ_units)

        no_outliers = remove_outliers(data)
        if length(λ) > 20
            smoothed, smoothed_err, _ = moving_average(no_outliers)
        else
            smoothed = no_outliers
            smoothed_err = zeros(size(smoothed))
        end

        if length(λ) > 20
            contrem, continuum = lunar_double_line(smoothed, λ)
        else
            contrem, continuum = single_line(smoothed, λ, (750, 1550))
        end

        spectral_resolution = (maximum(λ) - minimum(λ)) / length(λ)

        new{T, V}(
            data, λ, spectral_units, λ_units, spectral_resolution, no_outliers,
            smoothed, smoothed_err, contrem, continuum
        )
    end
end

Spectrum(
    data::AbstractVector{<:AbstractFloat},
    λ::AbstractVector{<:AbstractFloat},
    spectral_units::String,
    λ_units::String
) = Spectrum{Float64}(
    Float64.(data),
    Float64.(λ),
    spectral_units,
    λ_units
)

function validate_spectrum(data, λ, spectral_units, λ_units)

    size(data) == size(λ) || error(
        "Spectral Data and Wavelength Data are not the same size."
    )

    !any(isnan.(data)) || error(
        "NaN values present in the spectrum."
    )

    !any(isnan.(λ)) || error(
        "NaN values present in the wavelength vector."
    )

    valid_spectral_units = [
        "Reflectance",
        "RFL",
        "Radiance",
        "RDN"
    ]
    spectral_units in valid_spectral_units || error(
        "Invalid spectral unit passed."
    )

    valid_λ_units = [
        "nanometers",
        "nm",
        "microns",
        "um"
    ]
    λ_units in valid_λ_units || error(
        "Invalid wavelength unit passed."
    )
end
