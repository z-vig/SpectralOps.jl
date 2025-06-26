# continuum_removal/double_line.jl

"
    lunar_double_line(
        spectrum::AbstractVector{T},
        λ::AbstractVector{T}
    )::Tuple{AbstractVector{T}, AbstractVector{T}} where T<:AbstractFloat

Following method presented in Henderson et al., 2023
First, a rough continuum is removed using fixed points at 700, 1550 and 2600 nm
Next, three points are chosen from the maxima of this spectrum at:
 + 650 - 1000 nm
 + 1350 - 1600 nm
 + 2000 - 2600 nm
Finally, with these endpoints, the final continuum is calculated from the rfl
values at these points on the original spectrum

# Arguments
- `spectrum::AbstractVector{T}`: Input spectrum data.
- `λ::AbstractVector{T}`: Wavelength values of the spectrum.

# Returns
- `cont1_rem::AbstractVector{T}`: Spectrum with continuum removed.
- `cont1_complete::AbstractVector{T}`: Continuum values.
"
function lunar_double_line(
    spectrum::AbstractVector{T},
    λ::AbstractVector{T}
)::Tuple{Vector{Float64}, Vector{Float64}} where T <: AbstractFloat
    #Getting initial continuum line
    cont1_band_indices = [findλ(λ,i)[1] for i ∈ [700,1550,2600]]

    cont1_wvls = [findλ(λ,i)[2] for i ∈ [700,1550,2600]]
    cont1_spectrum_values = spectrum[cont1_band_indices]

    cont1_complete = linear_interpolation(
        cont1_wvls,
        cont1_spectrum_values,
        λ
    )

    cont1_rem = spectrum./cont1_complete

    # return cont1_complete,cont1_rem
    RANGE1 = (650,1000)
    RANGE2 = (1350,1600)
    RANGE3 = (2000,2600)

    cont2_band_indices = zeros(Int,3)
    n = 1
    for (i,j) ∈ [RANGE1,RANGE2,RANGE3]
        min_index = findλ(λ,i)[1]
        max_index = findλ(λ,j)[1]
        cont2_band_indices[n] = argmax(
            cont1_rem[range(min_index,max_index)]
        )+(min_index-1)
        n+=1
    end
    @debug cont2_band_indices

    cont2_wvls = λ[cont2_band_indices]
    cont2_spectrum_values = spectrum[cont2_band_indices]

    cont2_complete = linear_interpolation(
        cont2_wvls,
        cont2_spectrum_values,
        λ
    )

    cont2_rem = spectrum./cont2_complete

    return cont2_rem, cont2_complete
end

function lunar_double_line(
    cube::AbstractArray{T, 3},
    λ::AbstractVector{T}
)::Tuple{Array{Float64, 3}, Array{Float64, 3}} where T <: AbstractFloat
    ax = axes(cube)
    contrem_res = map(CartesianIndices(ax[1:2])) do I
        x,y = Tuple(I)
        spec = cube[x,y,:]
        return lunar_double_line(spec, λ)
    end

    contrem, continuum = separate_tuple_matrix(contrem_res)

    return contrem, continuum
end