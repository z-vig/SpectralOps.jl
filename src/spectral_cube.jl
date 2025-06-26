# Spectrum.jl

using HDF5

"
    SpectralCube(
        data::AbstractArray{T, 3}
        λ::AbstractVector{T}
        spectral_units::String
        λ_units::String
    )

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
struct SpectralCube{
    T<:AbstractFloat,
    A<:AbstractArray{T, 3},
    V<:AbstractVector{T}
}
    data::A
    λ::V
    spectral_units::String
    λ_units::String
    no_outliers::A
    smoothed::A
    smoothed_err::A
    contrem::A
    continuum::A

    function SpectralCube(
        data::A,
        λ::V,
        spectral_units::String,
        λ_units::String
    ) where {
        T<:AbstractFloat,
        A<:AbstractArray{T, 3},
        V<:AbstractVector{T}
    }
        validate_spectralcube(data, λ, spectral_units, λ_units)

        no_outliers = remove_outliers(data)
        println("Outliers removed.")

        if length(λ) < 20
            smoothed = no_outliers
            smoothed_err = zeros(size(no_outliers))
        else
            smoothed, smoothed_err, _ = moving_average(no_outliers)
        end
        println("Spectra smoothed.")

        if length(λ) < 20
            contrem, continuum = single_line(smoothed, λ, (750, 1550))
        else
            contrem, continuum = lunar_double_line(smoothed, λ)
        end

        println("Continuum removed.")

        new{T, A, V}(
            data, λ, spectral_units, λ_units, no_outliers, smoothed,
            smoothed_err, contrem, continuum
        )
    end

    function SpectralCube(
        data::A,
        λ::V,
        spectral_units::String,
        λ_units::String,
        no_outliers::A,
        smoothed::A,
        smoothed_err::A,
        contrem::A,
        continuum::A
    ) where {
        T<:AbstractFloat,
        A<:AbstractArray{T, 3},
        V<:AbstractVector{T}
    }
        new{T, A, V}(
            data, λ, spectral_units, λ_units, no_outliers, smoothed,
            smoothed_err, contrem, continuum
        )
    end
end

SpectralCube(
    data::AbstractArray{<:AbstractFloat, 3},
    λ::AbstractVector{<:AbstractFloat},
    spectral_units::String,
    λ_units::String
) = SpectralCube{Float64}(
    Float64.(data),
    Float64.(λ),
    spectral_units,
    λ_units
)

function validate_spectralcube(data, λ, spectral_units, λ_units)

    size(data)[end] == size(λ)[end] || error(
        "Spectral Data and Wavelength Data are not the same size."
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


save_data = Dict(
    :data => "original",
    :no_outliers => "outliers_removed",
    :smoothed => "smoothed",
    :smoothed_err => "smoothed_error",
    :contrem => "continuum_removed",
    :continuum => "continuum"
)

save_attrs = Dict(
    :λ => "wavelengths",
    :spectral_units => "spectral_units",
    :λ_units => "wavelength_units"
)

"
    save_spectral_cube(cube::SpectralCube, h5file::HDF5.File)

Saves SpectralCube object to an HDF5 file.
"
function save_spectral_cube(cube::SpectralCube, h5file::HDF5.File)
    for (key, val) ∈ save_data
        h5file[val] = getfield(cube, key)
    end

    for (key, val) ∈ save_attrs
        attrs(h5file)[val] = getfield(cube, key)
    end
end


"
    read_spectral_cube(h5path::String)::SpectralCube

Reads SpectralCube object from a saved HDF5 file.
"
function read_spectral_cube(h5path::String)::SpectralCube
    field_vals = []

    h5open(h5path, "r") do file
        for val ∈ values(save_attrs)
            push!(field_vals, attrs(file)[val])
        end

        for val ∈ values(save_data)
            push!(field_vals, read(file, val))
        end

        return nothing
    end

    field_vals = [field_vals[4], field_vals[1:3]..., field_vals[5:end]...]

    return SpectralCube(field_vals...)
end