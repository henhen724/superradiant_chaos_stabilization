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
    u[2:end] /= norm(u[2:end])
    ω_c = ω_tilde * E_0
    trans_sum_all = 0
    checker_board_all = 0
    for n in -transmax:transmax
        for m in -longmax:longmax
            mom_indx = to_1d_index(n, m, transmax, longmax)
            long_sum = safe_index_2D(u, n + 2, m, longmax, transmax) + safe_index_2D(u, n - 2, m, longmax, transmax)
            trans_sum = safe_index_2D(u, n, m + 2, longmax, transmax) + safe_index_2D(u, n, m - 2, longmax, transmax)
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
                H[1+to_1d_index(n, m, transmax, longmax), 1+to_1d_index(n + 1, m + 1, transmax, longmax)] = -ω_r * (P * conj(λ) + conj(P) * λ)
            end
            if n + 1 <= transmax && m - 1 >= -longmax
                H[1+to_1d_index(n, m, transmax, longmax), 1+to_1d_index(n + 1, m - 1, transmax, longmax)] = -ω_r * (P * conj(λ) + conj(P) * λ)
            end
            if n - 1 >= -transmax && m + 1 <= longmax
                H[1+to_1d_index(n, m, transmax, longmax), 1+to_1d_index(n - 1, m + 1, transmax, longmax)] = -ω_r * (P * conj(λ) + conj(P) * λ)
            end
            if n - 1 >= -transmax && m - 1 >= -longmax
                H[1+to_1d_index(n, m, transmax, longmax), 1+to_1d_index(n - 1, m - 1, transmax, longmax)] = -ω_r * (P * conj(λ) + conj(P) * λ)
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
            trans_sum = safe_index_2D(u, n, m + 2, longmax, transmax) + safe_index_2D(u, n, m - 2, longmax, transmax)
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
        u0 /= norm(u0)
    end

    function f(dx, x, p)
        λ = real_to_complex(x)
        dλ = real_to_complex(dx)
        dλ[1] = cavity_eq_for_eigvec(eigen_index, λ[1], u0; P=P, ω_tilde=ω_tilde, N_A=N_A, κ=κ, E_0=E_0, ω_r=ω_r, longmax=longmax, transmax=longmax)
        complex_to_real(dx, dλ)
    end

    prob = NonlinearProblem(f, complex_to_real([λ0]), nothing; abstol=1e-6, reltol=1e-6)
    sol = solve(prob, RobustMultiNewton(; autodiff=AutoFiniteDiff()); abstol=1e-6)

    @assert sol.retcode == ReturnCode.Success

    λ = real_to_complex(sol.u)[1]

    H = atomic_hamiltonian!(λ; P=P, ω_r=ω_r, longmax=longmax, transmax=transmax)

    eigvals, eigvecs = eigsolve(H, u0, eigen_index, :SR)
    vec = Array(eigvecs[eigen_index])

    return hcat([λ, vec...])[:, 1]
end

function find_steady_state_and_stability(; P=0.5, ω_tilde=1.5, N_A=10^5, κ=8.1, E_0=40.0, ω_r=0.05, longmax=2, transmax=2, eigen_index=1, λ0=3.0 + 3.0im)
    u_end = find_steady_state(; P=P, ω_tilde, N_A=N_A, κ=κ, ω_r=ω_r, longmax=longmax, transmax=transmax, eigen_index=eigen_index, λ0=λ0)

    function jacobian_multimomenta_model_drift!(J, u, p, t)
        J .= ForwardDiff.jacobian(temp_u -> begin
                uComplex = real_to_complex(temp_u)
                du = similar(uComplex)
                uComplex[2:end] /= norm(uComplex[2:end])
                multimomenta_model_drift!(du, uComplex, p, t; P=P, ω_tilde=ω_tilde, N_A=N_A, κ=κ, E_0=E_0, ω_r=ω_r, longmax=longmax, transmax=transmax)
                du[2:end] -= (du[2:end]' * uComplex[2:end]) / (uComplex[2:end]' * uComplex[2:end]) * uComplex[2:end]
                return complex_to_real(du)
            end, u)
    end

    # Calculate Jacobian at u_end
    A = zeros(Float64, 2 * length(u_end), 2 * length(u_end))
    jacobian_multimomenta_model_drift!(A, complex_to_real(u_end), nothing, 0.0)
    eigenvalues = eigvals(A)
    println(eigenvalues[real(eigenvalues).>=10^(-6)])
    is_stable = all(real(eigenvalues) .< 10^(-5))
    return u_end, is_stable, maximum(real.(eigenvalues))
end

function check_tol(a::T, tol) where {T<:Number}
    if abs(a) > tol
        return a
    end
    return zero(T)
end

function find_image_projection(M; tol=1e-6)
    U, s, Vt = svd(M)
    non_zero_indeces = [i for i in 1:size(U, 2) if abs(s[i]) > tol]
    return U[:, non_zero_indeces]
end

function gram_schmidt(V; tol=1e-6, check_from=2)
    U = copy(V)
    U[:, 1] /= norm(U[:, 1])
    for i in check_from:size(V, 2)
        U[:, i] -= U[:, 1:i-1] * (U[:, 1:i-1]' * U[:, i])
        if norm(U[:, i]) > tol
            U[:, i] /= norm(U[:, i])
        else
            U[:, i] .= 0.0
        end
    end
    return U
end

function remove_zero_columns!(matrix; tol=1e-6)
    non_zero_cols = [i for i in 1:size(matrix, 2) if norm(matrix[:, i]) > tol]
    return matrix[:, non_zero_cols]
end

function controllable_subspace_projector(A, B; n=size(A, 2), tol=1e-6)
    Pt = find_image_projection(Array(B); tol=tol)
    B_curr = B
    for i in 1:(n-1)
        B_curr = A * B_curr
        Pt = find_image_projection(Array(hcat(Pt, B_curr)); tol=tol)
        if rank(Pt) == size(A, 1)
            break
        end
    end
    return Pt
end

function controllable(A, B; tol=1e-10)
    Pt = controllable_subspace_projector(A, B; tol=tol)
    return rank(Pt) == size(A, 1)
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

function plot_complex_function(f; xlims=(-1.0, 1.0), ylims=(-1.0, 1.0), resolution=100, gamma=1, kwargs...)
    x = range(xlims[1], xlims[2], length=resolution)
    y = range(ylims[1], ylims[2], length=resolution)
    z = [f(complex(re, im)) for re in x, im in y]

    plot_complex_mesh(x, y, z; xlims=xlims, ylims=ylims, gamma=gamma, color=:hsv, kwargs...)
end

function plot_complex_mesh(x, y, z; gamma=0.5, kwargs...)
    angles = angle.(z)
    mags = (abs.(z)) .^ gamma
    colors = HSV.((angles) ./ (π / 180.0), 1.0, mags ./ maximum(mags))

    heatmap(x, y, colors; kwargs...)#, kwargs...)
end

function plot_gaussian_bumps(u0; longmax=2, transmax=2, sigma=0.1, pixel_per_bump=10, kwargs...)
    x = range(-(transmax + 1), (transmax + 1), length=pixel_per_bump * (2 * transmax + 1))
    y = range(-(longmax + 1), (longmax + 1), length=pixel_per_bump * (2 * longmax + 1))
    X, Y = [x for x in x, y in y], [y for x in x, y in y]
    Z = zeros(size(X))

    Threads.@threads for n in -transmax:transmax
        for m in -longmax:longmax
            idx = to_1d_index(n, m, transmax, longmax)
            amplitude = abs(u0[2+idx])
            Z .+= amplitude * exp.(-((X .- n) .^ 2 .+ (Y .- m) .^ 2) / (2 * sigma^2))
        end
    end

    plot_complex_mesh(x, y, Z; color=:viridis, xlabel=L"$\frac{k_x}{k_r}$", ylabel=L"$\frac{k_z}{k_r}$", xlims=(x[begin], x[end]), ylim=(y[begin], y[end]), kwargs...)
end