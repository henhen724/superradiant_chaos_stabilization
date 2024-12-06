using DifferentialEquations, Plots, LaTeXStrings, LinearAlgebra, KrylovKit, NonlinearSolve
include("../multimomenta_lib.jl")

P = 0.5
ω_tilde = 1.5
N_A = 10^5
κ = 8.1
E_0 = 40.0
ω_r = 0.05
longmax = 10
transmax = 10

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

λ = find_lambda(; λ0=2.0, P=P, ω_tilde=ω_tilde, N_A=N_A, κ=κ, E_0=E_0, ω_r=ω_r, longmax=longmax, transmax=transmax, n=1)