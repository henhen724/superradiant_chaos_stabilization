using DifferentialEquations, Plots, LaTeXStrings, NonlinearSolve, LinearAlgebra, NLsolve
include("multimomenta_lib.jl")

P = 0.0
ω_tilde = 3.00
N_A = 10^5
κ = 8.1
E_0 = 40.0
ω_r = 0.05
longmax = 10
transmax = 10

ω_c = E_0 * ω_tilde
# phi_e = sqrt(1 / 2 - ω_r / (4 * ω_c * P) * (ω_c^2 + κ^2))
phi_e = 0.4

vec_dim = 1 + (2 * longmax + 1) * (2 * transmax + 1)
u0 = (1 - im) / (2 * sqrt(2 * N_A)) * ones(ComplexF64, vec_dim)
u0[1] = 0.0
u0[2+to_1d_index(0, 0, transmax, longmax)] = sqrt(1 - phi_e^2)
u0[2+to_1d_index(1, 0, transmax, longmax)] = phi_e / 2
u0[2+to_1d_index(0, 1, transmax, longmax)] = phi_e / 2
u0[2+to_1d_index(-1, 0, transmax, longmax)] = phi_e / 2
u0[2+to_1d_index(0, -1, transmax, longmax)] = phi_e / 2
u0norm = sum(abs.(u0[2:end]) .^ 2)
u0[2:end] = u0[2:end] / sqrt(u0norm)


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

duComplex = 0.0 * similar(u0)

function fixed_point_condition(du, u, p, t)
    u[2:end] = u[2:end] / norm(u[2:end])
    # print(norm(u[2:end]))
    multimomenta_model_drift!(duComplex, u, nothing, nothing; P=P, ω_tilde=ω_tilde, N_A=N_A, κ=κ, E_0=E_0, ω_r=ω_r, longmax=longmax, transmax=transmax)

    du[1] = duComplex[1]

    du[2:end] = duComplex[2:end] + (1 - dot(conj.(u[2:end]), u[2:end]) - dot(conj.(duComplex[2:end]), u[2:end]) / dot(conj.(u[2:end]), u[2:end])) * u[2:end]

    return du
end

# u0 = sol.u

# prob = SteadyStateProblem(fixed_point_condition, u0)

# sol = solve(prob; abstol=1e-5, reltol=1e-5, maxiters=10^10, progress_steps=true)

