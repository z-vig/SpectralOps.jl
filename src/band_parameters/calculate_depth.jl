# band_parameters/depth.jl

"
    calculate_depth(
        spectrum::AbstractVector{T},
        λmin::T,
        λmax::T
    )::T where T <: AbstractFloat

Calculates the depth of a spectral absorption band.

# Arguments
- `spectrum::AbstractVector{T}`: Spectral data containing an absorption.
- `λmin::T`: Minimum wavelength (matching units of spectrum) to use in
             calculation.
- `λmax::T`: Maximum wavelength to use in calculation.

# Returns
- `depth::Float64`: Absorption depth at the absorption center.
"
function calculate_depth(
    contrem_spectrum::AbstractVector{T},
    λ::AbstractVector{T},
    λmin::V,
    λmax::V
)::Float64 where {
    T<:AbstractFloat,
    V<:Real
}
    center, xfit, yfit = calculate_center(contrem_spectrum, λ, λmin, λmax)
    depth = 1 - yfit[findλ(xfit, center)[1]]
    if depth < 1
        return NaN64
    else
        return depth
    end
end


function calculate_depth(
    contrem_cube::AbstractArray{T},
    λ::AbstractVector{T},
    λmin::V,
    λmax::V
)::Matrix{Float64} where {
    T<:AbstractFloat,
    V<:Real
}
    bc_map, xfit, yfit_cube = calculate_center(contrem_cube, λ, λmin, λmax)
    depth_map = Array{Float64}(undef, size(bc_map)...)
    for i ∈ axes(bc_map, 1)
        for j ∈ axes(bc_map, 2)
            depth_map[i, j] = 1 - view(
                yfit_cube, i, j, :
            )[findλ(xfit, bc_map[i, j])[1]]
        end
    end
    depth_map[depth_map .< 0] .= NaN64
    depth_map[depth_map .> 0.5] .= NaN64
    return depth_map
end


function calculate_depth(
    spec::Spectrum{T},
    λmin::V,
    λmax::V
)::Float64 where {
    T<:AbstractFloat,
    V<:Real
}
    calculate_depth(
        spec.contrem,
        spec.λ,
        λmin,
        λmax
    )
end

function calculate_depth(
    spec::SpectralCube{T},
    λmin::V,
    λmax::V
)::Matrix{Float64} where {
    T<:AbstractFloat,
    V<:Real
}
    calculate_depth(
        spec.contrem,
        spec.λ,
        λmin,
        λmax
    )
end
