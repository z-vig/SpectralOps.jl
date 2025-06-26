"
    seperate_tuple_matrix(
        matrix::Matrix{<:Tuple{<:Vector{T1}, <:Vector{T1}, <:Vector{T2}}}
    )::Tuple{Matrix} where {T1<:Real, T2<:Real}

Transforms a matrix of vector tuples into a tuple of seperate 3D matrices.
"
function separate_tuple_matrix(
    matrix::Matrix{
        <:Tuple{Vararg{Vector}}
    }
)::Tuple{Vararg{Array}}
    # Determine the size of the input matrix
    rows, cols = size(matrix)
    
    # Number of Vectors in each Tuple in the Matrix
    N = length(matrix[1, 1])

    # Initialize N empty matrices to hold the separated vectors
    separated_matrices = [
        Matrix{Vector{Float64}}(undef, rows, cols) for _ ∈ 1:N
    ]
    
    # Iterate over the input matrix and populate the new matrices
    for i in 1:rows
        for j in 1:cols
            vectuple = matrix[i, j]
            @assert length(vectuple) == N "All tuples must have the same length"
            for k ∈ 1:N
                separated_matrices[k][i, j] = vectuple[k]
            end
        end
    end
    
    separated_matrices_3D = [
        make3d(mat) for mat ∈ separated_matrices
    ]

    return tuple(separated_matrices_3D...)
end