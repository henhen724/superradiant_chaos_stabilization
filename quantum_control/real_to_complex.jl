function complex_to_real(vec::Vector{Complex{T}}) where {T}
    vec_dim = length(vec)
    vecReal = zeros(T, 2 * vec_dim)
    vecReal[begin:vec_dim] = real.(vec)
    vecReal[vec_dim+1:end] = imag.(vec)
    return vecReal
end

function complex_to_real(vecReal::Vector{T}, vec::Vector{Complex{T}}) where {T}
    vec_dim = length(vec)
    vecReal[begin:vec_dim] = real.(vec)
    vecReal[vec_dim+1:end] = imag.(vec)
end

function real_to_complex(vec::Vector{T}) where {T<:Real}
    @assert length(vec) % 2 == 0
    vec_dim = Int(length(vec) // 2)
    vecComplex = zeros(Complex{T}, vec_dim)
    vecComplex = complex.(vec[begin:vec_dim], vec[vec_dim+1:end])
    return vecComplex
end

function real_to_complex(vecComplex::Vector{Complex{T}}, vec::Vector{T}) where {T<:Real}
    @assert length(vec) % 2 == 0
    vec_dim = Int(length(vec) // 2)
    vecComplex .= complex.(vec[begin:vec_dim], vec[vec_dim+1:end])
end

function real_to_complex(mat::Matrix{T}) where {T<:Real}
    dims = size(mat)
    @assert dims[1] % 2 == 0 && dims[2] % 2 == 0
    c_dims = (Int(dims[1] // 2), Int(dims[2] // 2))
    c_matrix = zeros(Complex{T}, (dims[1], dims[2]))
    c_matrix[begin:c_dims[1], begin:c_dims[2]] = complex.(mat[begin:c_dims[1], begin:c_dims[2]] + mat[c_dims[1]+1:end, c_dims[2]+1:end], mat[c_dims[1]+1:end, begin:c_dims[2]] - mat[begin:c_dims[1], c_dims[2]+1:end])
    c_matrix[c_dims[1]+1:end, begin:c_dims[2]] = complex.(mat[begin:c_dims[1], begin:c_dims[2]] + mat[c_dims[1]+1:end, c_dims[2]+1:end], -mat[c_dims[1]+1:end, begin:c_dims[2]] + mat[begin:c_dims[1], c_dims[2]+1:end])
    c_matrix[begin:c_dims[1], c_dims[2]+1:end] = complex.(mat[begin:c_dims[1], begin:c_dims[2]] - mat[c_dims[1]+1:end, c_dims[2]+1:end], mat[c_dims[1]+1:end, begin:c_dims[2]] + mat[begin:c_dims[1], c_dims[2]+1:end])
    c_matrix[c_dims[1]+1:end, c_dims[2]+1:end] = complex.(mat[begin:c_dims[1], begin:c_dims[2]] - mat[c_dims[1]+1:end, c_dims[2]+1:end], -mat[c_dims[1]+1:end, begin:c_dims[2]] - mat[begin:c_dims[1], c_dims[2]+1:end])
    return c_matrix ./ 2
end