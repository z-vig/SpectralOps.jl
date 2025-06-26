# utils/polyfit.jl

"
    vandermonde(
        x::AbstractVector{T},
        order::Int
    )::Matrix{Float64} where T<:AbstractFloat

Returns vandermonde matrix of `order` N using `x` data.
"
function vandermonde(
    x::AbstractVector{T},
    order::Int
)::Matrix{Float64} where T<:AbstractFloat
    return reduce(hcat, (x.^i for i in 0:order))
end

"
    polyfit(
        x::Union{Nothing, Vector{T}},
        y::Vector{T},
        order::Union{Nothing, Int},
        design_matrix::Matrix{T}
    )::Vector{Float64}

Fits a polynomial of `order` N to `x` and `y` data. Optionally can be used in a
loop by specifying the `design_matrix` elements.

# Arguments
- `x::Union{Nothing, Vector{T}}`: X Data.
- `y::AbstractVector{T}`: Y Data.
- `order::Union{Nothing, Int}`: Order of polynomial.
- `design_matrices::Union{Nothing, NTuple{3, Matrix{T}}`: Optional. Precomputed
        design matrix components, (X, Xt and XtX). If None (default), these
        are calculated from `x` data.

# Results
- `β::Vector{Float64}`: Coefficients of the fit.
"
function polyfit(
    x::Union{Nothing, AbstractVector{T}},
    y::AbstractVector{T},
    order::Union{Nothing, Int};
    design_matrices::Union{Nothing, NTuple{3, AbstractMatrix{T}}} = nothing
)::Vector{Float64} where T<:AbstractFloat
    if isnothing(design_matrices)
        if isnothing(x) || isnothing(order)
            error("If desigin matrices are not given, x data and order must be
                  given.")
        end
        X = vandermonde(x, order)
        XtX = X' * X
        Xt = X'
    else
        X, Xt, XtX = design_matrices
    end


    β = XtX \ (Xt * y)

    return β
end

"
    polyfit(
        spectral_cube::AbstractArray{T, 3},
        λ::AbstractVector{T},
        order::Int
    )::Matrix{Float64}

Fits polynomials to an entire cube of spectra.

# Arguments
- `spectral_cube::AbstractArray{T, 3}`: Spectral image cube.
- `λ::AbstractVector{T}`: Wavelength vector.
- `order::Int`: Order of polynomial fit.
"
function polyfit(
    spectral_cube::AbstractArray{T, 3},
    λ::AbstractVector{T},
    order::Int
)::Matrix{Float64} where T<:AbstractFloat
    X = vandermonde(λ, order)
    Xt = X'
    XtX = Xt * X
    design_matrices = (X, Xt, XtX)
    
    N, M, B = size(spectral_cube)
    nspectra = N*M
    spectra_to_fit = reshape(spectral_cube, nspectra, B)

    coefs = Array{Float64}(undef, nspectra, size(X, 2))

    for i ∈ 1:nspectra
        y = @view spectra_to_fit[i, :]
        coefs[i, :] = polyfit(nothing, y, nothing; design_matrices)
    end
    
    return coefs
end
