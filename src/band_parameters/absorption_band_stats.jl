# band_parameters/absorption_band_stats.jl

"
    AbsorptionBandStats(
        spectrum::AbstractVector{T},
        λ_min::T,
        λ_max::T
    )

Calculates and stores information about a single absorption band within a
single spectrum.
"
struct AbsorptionBandStats{T, V}
    spectrum::Spectrum{T}
    λ_min::V
    λ_max::V
    depth::T
    center::T
    area::T
    asymmetry::T

    function AbsorptionBandStats(
        spectrum::Spectrum{T},
        λ_min::V,
        λ_max::V
    ) where {
        T<:AbstractFloat,
        V<:Real
    }
        area = calculate_area(spectrum, λ_min, λ_max)
        depth = calculate_depth(spectrum, λ_min, λ_max)
        center = calculate_center(spectrum.data, λ_min, λ_max)
        asymmetry = calculate_asymmetry(spectrum.data, λ_min, λ_max)

        new{T, V}(spectrum, λ_min, λ_max, area, depth, center, asymmetry)
    end
end
