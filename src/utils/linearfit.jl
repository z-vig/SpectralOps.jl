# utils/linearfit.jl

"
    linearfit(
        x::AbstractVector{T},
        y::AbstractVector{T},
        xfit::AbstractVector{T}
    ):: T where T<:AbstractFloat

Performs a linear fit to data.

# Arguments
- `x::AbstractVector{T}`: X-data.
- `y::AbstractVector{T}`: Y-data.
- `xfit::AbstractVector{T}`: X data to use in fitted line.

# Returns
- `yfit::AbstractVector{T}`: Fitted y-data.
"
function linearfit(
    x::AbstractVector{T1},
    y::AbstractVector{T2},
    xfit::AbstractVector{T3}
):: Vector{Float64} where {T1<:Real, T2<:Real, T3<:Real}

    # Building design matrix
    G = Array{Float64}(undef, length(x), 2)
    G[:, 1] .= 1
    G[:, 2] = x

    # Building data matrix
    d = Array{Float64}(undef, length(y), 1)
    d[:, 1] = y

    m = G' * G
    m = inv(m)
    m = m * G'
    m = m*d

    intercept = m[1, 1]
    slope = m[2, 1]
    yfit = slope .* xfit .+ intercept

    return yfit
end
