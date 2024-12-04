using SparseArrays

function to_1d_index(n::Int, m::Int, transmax::Int, longmax::Int)::Int
    @assert -transmax <= n <= transmax
    @assert -longmax <= m <= longmax
    return (n + transmax) * (2 * longmax + 1) + m + longmax
end

function safe_index_2D(vec::Vector{T}, n::Int, m::Int, transmax::Int, longmax::Int) where {T}
    if -transmax <= n <= transmax && -longmax <= m <= longmax
        return vec[2+to_1d_index(n, m, transmax, longmax)]
    else
        return zero(T)
    end
end

function multimomenta_model_drift!(du, u, p, t; P=0.0, ω_tilde=0.0, N_A=10^5, κ=8.1, E_0=40.0, ω_r=0.05, longmax=10, transmax=10)
    ω_c = ω_tilde * E_0
    trans_sum_all = 0
    checker_board_all = 0
    for n in -transmax:transmax
        for m in -longmax:longmax
            mom_indx = to_1d_index(n, m, transmax, longmax)
            trans_sum = safe_index_2D(u, n + 2, m, longmax, transmax) + safe_index_2D(u, n - 2, m, longmax, transmax)
            long_sum = safe_index_2D(u, n, m + 2, longmax, transmax) + safe_index_2D(u, n, m - 2, longmax, transmax)
            trans_sum_all += conj(trans_sum) * u[2+mom_indx]
            checker_board = safe_index_2D(u, n + 1, m + 1, longmax, transmax) +
                            safe_index_2D(u, n + 1, m - 1, longmax, transmax) +
                            safe_index_2D(u, n - 1, m + 1, longmax, transmax) +
                            safe_index_2D(u, n - 1, m - 1, longmax, transmax)
            checker_board_all += conj(checker_board) * u[2+mom_indx]
            du[2+mom_indx] = -im * ω_r * ((n^2 + m^2) * u[2+mom_indx] - conj(u[1]) * u[1] * long_sum - P * (u[1] + conj(u[1])) * checker_board - P^2 * trans_sum)
        end
    end
    du[1] = -(κ + im * ω_c) * u[1] + im * E_0 * u[1] * trans_sum_all + im * E_0 * P * checker_board_all
end

function atomic_hamiltonian!(λ; P=0.0, ω_r=0.05, longmax=10, transmax=10)
    vec_dim = (2 * longmax + 1) * (2 * transmax + 1)
    H = spzeros(ComplexF64, (vec_dim, vec_dim))
    for n in -transmax:transmax
        for m in -longmax:longmax
            H[1+to_1d_index(n, m, transmax, longmax), 1+to_1d_index(n, m, transmax, longmax)] = ω_r
            if n + 2 <= transmax
                H[1+to_1d_index(n, m, transmax, longmax), 1+to_1d_index(n + 2, m, transmax, longmax)] = -ω_r * conj(λ) * λ
            end
            if n - 2 >= -transmax
                H[1+to_1d_index(n, m, transmax, longmax), 1+to_1d_index(n - 2, m, transmax, longmax)] = -ω_r * conj(λ) * λ
            end
            if m + 2 <= longmax
                H[1+to_1d_index(n, m, transmax, longmax), 1+to_1d_index(n, m + 2, transmax, longmax)] = -ω_r * P^2
            end
            if m - 2 >= -longmax
                H[1+to_1d_index(n, m, transmax, longmax), 1+to_1d_index(n, m - 2, transmax, longmax)] = -ω_r * P^2
            end
            if n + 1 <= transmax && m + 1 <= longmax
                H[1+to_1d_index(n, m, transmax, longmax), 1+to_1d_index(n + 1, m + 1, transmax, longmax)] = -ω_r * P * (conj(λ) + λ)
            end
            if n + 1 <= transmax && m - 1 >= -longmax
                H[1+to_1d_index(n, m, transmax, longmax), 1+to_1d_index(n + 1, m - 1, transmax, longmax)] = -ω_r * P * (conj(λ) + λ)
            end
            if n - 1 >= -transmax && m + 1 <= longmax
                H[1+to_1d_index(n, m, transmax, longmax), 1+to_1d_index(n - 1, m + 1, transmax, longmax)] = -ω_r * P * (conj(λ) + λ)
            end
            if n - 1 >= -transmax && m - 1 >= -longmax
                H[1+to_1d_index(n, m, transmax, longmax), 1+to_1d_index(n - 1, m - 1, transmax, longmax)] = -ω_r * P * (conj(λ) + λ)
            end
        end
    end
    return dropzeros(H)
end

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
    vecComplex = vec[begin:vec_dim] + im * vec[vec_dim+1:end]
    return vecComplex
end

function real_to_complex(vecComplex::Vector{Complex{T}}, vec::Vector{T}) where {T<:Real}
    @assert length(vec) % 2 == 0
    vec_dim = Int(length(vec) // 2)
    vecComplex .= vec[begin:vec_dim] + im * vec[vec_dim+1:end]
end