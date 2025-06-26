# utils/findwvl.jl

"""
    findλ(λ, targetλ)

Given a list of wavelengths, `λ`, find the index of a `targetλ` and the actual
wavelength closest to your target.
"""
function findλ(
    λ::Vector{T},
    targetλ::Real
)::Tuple{Int,Float64} where {T<:AbstractFloat}
    idx = argmin(abs.(λ .- targetλ))
    return (idx,λ[idx])
end