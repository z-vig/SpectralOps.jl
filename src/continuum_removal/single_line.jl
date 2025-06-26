"
    single_line(
        spectrum::AbstractVector{T},
        λ::AbstractVector{T}
    )::Tuple{AbstractVector{T}, AbstractVector{T}} where T<:AbstractFloat

Removes a single line from the spectrum as an estimate for the continuum.

# Arguments
- `spectrum::AbstractVector{T}`: Input spectrum data.
- `λ::AbstractVector{T}`: Wavelength values of the spectrum.
- `tie_bands::Tuple{Vararg{R}}`: Wavelengths to tie continuum to.

# Returns
- `contrem::AbstractVector{T}`: Spectrum with continuum removed.
- `continuum::AbstractVector{T}`: Continuum values.
"
function single_line(
    spectrum::AbstractVector{T},
    λ::AbstractVector{T},
    tie_bands::Tuple{Vararg{R}}
)::Tuple{Vector{Float64}, Vector{Float64}} where {
    T<:AbstractFloat,
    R<:Real
}
    #Getting initial continuum line
    continuum_indices = [findλ(λ,i)[1] for i ∈ tie_bands]

    continuum_wvls = [findλ(λ,i)[2] for i ∈ tie_bands]
    continuum_spectrum_values = spectrum[continuum_indices]

    continuum = linear_interpolation(
        continuum_wvls,
        continuum_spectrum_values,
        λ
    )

    contrem = spectrum./continuum

    return contrem, continuum
end

function single_line(
    cube::AbstractArray{T},
    λ::AbstractVector{T},
    tie_bands::Tuple{Vararg{R}}
)::Tuple{Array{Float64}, Array{Float64}} where {
    T<:AbstractFloat,
    R<:Real
}
    ax = axes(cube)
    contrem_res = map(CartesianIndices(ax[1:2])) do I
        x,y = Tuple(I)
        spec = cube[x,y,:]
        return single_line(spec, λ, tie_bands)
    end
    # Separate the matrix of tuples into three matrices
    contrem, continuum = separate_tuple_matrix(contrem_res)

    return contrem, continuum
end