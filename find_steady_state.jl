using DifferentialEquations, Plots, LaTeXStrings, LinearAlgebra, NLsolve, KrylovKit
include("multimomenta_lib.jl")

P = 0.5
ω_tilde = 1.5
N_A = 10^5
κ = 8.1
E_0 = 40.0
ω_r = 0.05
longmax = 10
transmax = 10

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

function find_lambda(; λ0=0.0, P=0.0, ω_tilde=0.0, N_A=10^5, κ=8.1, E_0=40.0, ω_r=0.05, longmax=10, transmax=10, n=1, rate=1.0)
    λ = λ0
    tol = 1e-5
    max_iters = 10^5
    iter = 0

    vec_dim = (2 * longmax + 1) * (2 * transmax + 1)
    u0 = (1 - im) / (2 * sqrt(2 * N_A)) * ones(ComplexF64, vec_dim)
    u0[1+to_1d_index(0, 0, transmax, longmax)] = 1.0
    u0[1+to_1d_index(1, 0, transmax, longmax)] = 0.0
    u0[1+to_1d_index(0, 1, transmax, longmax)] = 0.0
    u0[1+to_1d_index(-1, 0, transmax, longmax)] = 0.0
    u0[1+to_1d_index(0, -1, transmax, longmax)] = 0.0
    u0norm = sum(abs.(u0[2:end]) .^ 2)
    u0 = u0 / sqrt(u0norm)

    while iter < max_iters
        iter += 1

        cavity_value = cavity_eq_for_eigvec(n, λ, u0; P=P, ω_tilde=ω_tilde, N_A=N_A, κ=κ, E_0=E_0, ω_r=ω_r, longmax=longmax, transmax=longmax)


        if abs(cavity_value) < tol
            break
        end

        cavity_value_per = cavity_eq_for_eigvec(n, λ + tol, u0; P=P, ω_tilde=ω_tilde, N_A=N_A, κ=κ, E_0=E_0, ω_r=ω_r, longmax=longmax, transmax=longmax)

        # Update λ using Newton-Raphson method
        derivative = (cavity_value_per - cavity_value) / tol
        delta_lambda = rate * cavity_value / derivative
        if abs(delta_lambda) > 1e-1
            λ -= rate * delta_lambda
        else
            println("Using final steps")
            λ -= delta_lambda
        end
        println("New lambda $(λ) with resid $(norm(cavity_value)) and derivative $(derivative)")
    end

    if iter == max_iters
        println("Warning reached max iters.")
    end

    return λ
end

λ = find_lambda(; λ0=2.0, P=P, ω_tilde=ω_tilde, N_A=N_A, κ=κ, E_0=E_0, ω_r=ω_r, longmax=longmax, transmax=transmax, n=1, rate=0.1)


function julians_method(; λ0=0.0, P=0.0, ω_tilde=0.0, N_A=10^5, κ=8.1, E_0=40.0, ω_r=0.05, longmax=10, transmax=10, n=1, dt=1.0)
    λ = λ0
    tol = 1e-5
    max_iters = 10^5
    iter = 0

    vec_dim = (2 * longmax + 1) * (2 * transmax + 1)
    u0 = (1 - im) / (2 * sqrt(2 * N_A)) * ones(ComplexF64, vec_dim)
    u0[1+to_1d_index(0, 0, transmax, longmax)] = 1.0
    u0[1+to_1d_index(1, 0, transmax, longmax)] = 0.0
    u0[1+to_1d_index(0, 1, transmax, longmax)] = 0.0
    u0[1+to_1d_index(-1, 0, transmax, longmax)] = 0.0
    u0[1+to_1d_index(0, -1, transmax, longmax)] = 0.0
    u0norm = sum(abs.(u0[2:end]) .^ 2)
    u0 = u0 / sqrt(u0norm)

    while iter < max_iters
        iter += 1

        cavity_value = cavity_eq_for_eigvec(n, λ, u0; P=P, ω_tilde=ω_tilde, N_A=N_A, κ=κ, E_0=E_0, ω_r=ω_r, longmax=longmax, transmax=longmax)


        if abs(cavity_value) < tol
            break
        end

        λ -= dt * cavity_value
        println("New lambda $(λ) with resid $(norm(cavity_value))")
    end

    if iter == max_iters
        println("Warning reached max iters.")
    end

    return λ
end

λ = julians_method(; λ0=2.0, P=P, ω_tilde=ω_tilde, N_A=N_A, κ=κ, E_0=E_0, ω_r=ω_r, longmax=longmax, transmax=transmax, n=1, dt=10^(-6))
