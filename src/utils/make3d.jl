# utils/make3d.jl

"""
    make3d(
        matrix::AbstractMatrix{AbstractVector{T}}
    )::Array{T, 3} where T<:AbstractFloat

A function for turning a Matrix{Vector} to an Array{,3}
"""
function make3d(matrix::AbstractMatrix{V})::Array{T, 3} where {
    T <: AbstractFloat,
    V <: AbstractVector{T}
}
    return permutedims(
        [
            matrix[I][k]
            for k = eachindex(matrix[1,1]), I=CartesianIndices(matrix)
        ],
        (2,3,1)
    )
end