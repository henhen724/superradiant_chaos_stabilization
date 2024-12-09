using SparseArrays, KrylovKit

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

function multimomenta_model_drift!(du, u, p, t; P=0.0, ω_tilde=0.0, N_A=10^5, κ=8.1, E_0=40.0, ω_r=0.05, ϵ=0.0, longmax=10, transmax=10)
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
            du[2+mom_indx] = -im * ω_r * ((n^2 + m^2) * u[2+mom_indx] - conj(u[1]) * u[1] * long_sum - (conj(P) * u[1] + P * conj(u[1])) * checker_board - conj(P) * P * trans_sum)
        end
    end
    du[1] = ϵ - (κ + im * ω_c) * u[1] + im * E_0 * u[1] * trans_sum_all + im * E_0 * P * checker_board_all
end

function dispative_dynamics!(du, u, p, t; P=0.0, ω_tilde=0.0, N_A=10^5, κ=8.1, E_0=40.0, ω_r=0.05, longmax=10, transmax=10)
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
            du[2+mom_indx] = -ω_r * ((n^2 + m^2) * u[2+mom_indx] - conj(u[1]) * u[1] * long_sum - P * (u[1] + conj(u[1])) * checker_board - P^2 * trans_sum)
        end
    end
    du[1] = -(κ + im * ω_c) * u[1] + im * E_0 * u[1] * trans_sum_all + im * E_0 * P * checker_board_all
end

function atomic_hamiltonian!(λ::T; P=0.0, ω_r=0.05, longmax=10, transmax=10) where {T}
    vec_dim = (2 * longmax + 1) * (2 * transmax + 1)
    H = spzeros(T, (vec_dim, vec_dim))
    for n in -transmax:transmax
        for m in -longmax:longmax
            H[1+to_1d_index(n, m, transmax, longmax), 1+to_1d_index(n, m, transmax, longmax)] = ω_r * (n^2 + m^2)
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

function cavity_eq!(u; P=0.0, ω_tilde=0.0, N_A=10^5, κ=8.1, E_0=40.0, ω_r=0.05, longmax=10, transmax=10)
    ω_c = ω_tilde * E_0
    trans_sum_all = 0
    checker_board_all = 0
    for n in -transmax:transmax
        for m in -longmax:longmax
            mom_indx = to_1d_index(n, m, transmax, longmax)
            trans_sum = safe_index_2D(u, n + 2, m, longmax, transmax) + safe_index_2D(u, n - 2, m, longmax, transmax)
            trans_sum_all += conj(trans_sum) * u[2+mom_indx]
            checker_board = safe_index_2D(u, n + 1, m + 1, longmax, transmax) +
                            safe_index_2D(u, n + 1, m - 1, longmax, transmax) +
                            safe_index_2D(u, n - 1, m + 1, longmax, transmax) +
                            safe_index_2D(u, n - 1, m - 1, longmax, transmax)
            checker_board_all += conj(checker_board) * u[2+mom_indx]
        end
    end
    return (-(κ + im * ω_c) * u[1] + im * E_0 * u[1] * trans_sum_all + im * E_0 * P * checker_board_all)
end

function cavity_eq_for_eigvec(n, λ, u0; P=0.0, ω_tilde=0.0, N_A=10^5, κ=8.1, E_0=40.0, ω_r=0.05, longmax=10, transmax=10)
    H = atomic_hamiltonian!(λ; P=P, ω_r=ω_r, longmax=longmax, transmax=transmax)

    # Find the lowest energy eigenvector using eigsolve
    eigenvalues, eigenvectors = eigsolve(H, u0, n, :SR)
    lowest_energy_eigenvector = Array(eigenvectors[n])

    # Normalize the eigenvector
    lowest_energy_eigenvector /= norm(lowest_energy_eigenvector)

    return cavity_eq!(cat(λ, lowest_energy_eigenvector, dims=1); P=P, ω_tilde=ω_tilde, N_A=N_A, κ=κ, E_0=E_0, ω_r=ω_r, longmax=longmax, transmax=transmax)
end

function find_steady_state(; P=0.0, ω_tilde=2.0, N_A=10^5, κ=8.1, E_0=40.0, ω_r=0.05, longmax=10, transmax=10, eigen_index=1, u0=nothing, λ0=1.0 + im)
    if u0 isa Nothing
        vec_dim = (2 * longmax + 1) * (2 * transmax + 1)
        u0 = (1 - im) / (2 * sqrt(2 * N_A)) * ones(ComplexF64, vec_dim)
        u0[1+to_1d_index(0, 0, transmax, longmax)] = 1.0
        u0[1+to_1d_index(1, 0, transmax, longmax)] = 0.0
        u0[1+to_1d_index(0, 1, transmax, longmax)] = 0.0
        u0[1+to_1d_index(-1, 0, transmax, longmax)] = 0.0
        u0[1+to_1d_index(0, -1, transmax, longmax)] = 0.0
        u0norm = sum(abs.(u0[2:end]) .^ 2)
        u0 = u0 / sqrt(u0norm)
    end

    function f(dx, x, p)
        λ = real_to_complex(x)
        dλ = real_to_complex(dx)
        dλ[1] = cavity_eq_for_eigvec(eigen_index, λ[1], u0; P=P, ω_tilde=ω_tilde, N_A=N_A, κ=κ, E_0=E_0, ω_r=ω_r, longmax=longmax, transmax=longmax)
        complex_to_real(dx, dλ)
    end

    prob = NonlinearProblem(f, complex_to_real([λ0]), nothing; abstol=1e-3, reltol=1e-3)
    sol = solve(prob, RobustMultiNewton(; autodiff=AutoFiniteDiff()); abstol=1e-5)

    @assert sol.retcode == ReturnCode.Success

    λ = real_to_complex(sol.u)[1]

    H = atomic_hamiltonian!(λ; P=P, ω_r=ω_r, longmax=longmax, transmax=transmax)

    eigvals, eigvecs = eigsolve(H, u0, eigen_index, :SR)
    vec = Array(eigvecs[eigen_index])

    return hcat([λ, vec...])[:, 1]
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