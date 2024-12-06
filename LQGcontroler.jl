using LinearAlgebra, DifferentialEquations, ForwardDiff
include("multimomenta_lib.jl")

# Define parameters
P = 0.5
ω_tilde = 1.5
N_A = 10^5
κ = 8.1
E_0 = 40.0
ω_r = 0.05
longmax = 10
transmax = 10

# Initial conditions
vec_dim = 1 + (2 * longmax + 1) * (2 * transmax + 1)
u0 = (1 - im) / (2 * sqrt(2 * N_A)) * ones(ComplexF64, vec_dim)
u0[1] = 3.0
u0[2+to_1d_index(0, 0, transmax, longmax)] = 1.0
u0[2+to_1d_index(1, 0, transmax, longmax)] = 0.0
u0[2+to_1d_index(0, 1, transmax, longmax)] = 0.0
u0[2+to_1d_index(-1, 0, transmax, longmax)] = 0.0
u0[2+to_1d_index(0, -1, transmax, longmax)] = 0.0
u0norm = sum(abs.(u0[2:end]) .^ 2)
u0[2:end] = u0[2:end] / sqrt(u0norm)

# Define drift function
function drift!(du, u, p, t)
    u[2:end] = u[2:end] / norm(u[2:end])
    dispative_dynamics!(du, u, p, t; P=P, ω_tilde=ω_tilde, N_A=N_A, κ=κ, E_0=E_0, ω_r=ω_r, longmax=longmax, transmax=transmax)
    du[2:end] -= u[2:end] * dot(conj.(u[2:end]), du[2:end]) / dot(conj.(u[2:end]), u[2:end])
end

# Time span
trecord = 0.0:0.05:5000.0
tspan = (trecord[begin], trecord[end])

# Solve ODE
probODE = ODEProblem(drift!, u0, tspan)
sol = solve(probODE, Tsit5(); reltol=10^-4, abstol=10^-4, dt=10^(-3), maxiters=10^12, save_noise=false, save_everystep=false, saveat=trecord)

# End state
u_end = sol[end]
# Define Jacobian function
function jacobian_multimomenta_model_drift!(J, u, p, t)
    J .= ForwardDiff.jacobian(temp_u -> begin
            uComplex = real_to_complex(temp_u)
            du = similar(uComplex)
            multimomenta_model_drift!(du, uComplex, p, t; P=P, ω_tilde=ω_tilde, N_A=N_A, κ=κ, E_0=E_0, ω_r=ω_r, longmax=longmax, transmax=transmax)
            return complex_to_real(du)
        end, u)
end

# Calculate Jacobian at u_end
A = zeros(Float64, 2 * length(u_end), 2 * length(u_end))
jacobian_multimomenta_model_drift!(A, complex_to_real(u_end), nothing, 0.0)

# Find eigenvalues
eigenvalues = eigvals(A)

# Determine stability
is_stable = all(real(eigenvalues) .< 0)

println("Eigenvalues: ", eigenvalues)
println("Is the fixed point stable? ", is_stable)

# Define control matrix function
function control_matrix!(B, st_control, p, t)
    B .= ForwardDiff.jacobian(control -> begin
            cont_complex = real_to_complex(control)
            du = similar(u_end, eltype(cont_complex))
            multimomenta_model_drift!(du, u_end, p, t; P=cont_complex[1], ω_tilde=cont_complex[2], N_A=N_A, κ=κ, E_0=E_0, ω_r=ω_r, ϵ=cont_complex[3], longmax=longmax, transmax=transmax)
            return complex_to_real(du)
        end, st_control)
end

control_vec = ComplexF64[P, ω_tilde, 0.0]
B = zeros(Float64, 2 * length(u_end), 2 * 3)
control_matrix!(B, complex_to_real(control_vec), nothing, nothing)
Q = Array{Float64}(I(2 * length(u_end)))
R = Array{Float64}(I(6))

function solve_riccati(A, B, Q, R)
    P = lyap(A, B * inv(R) * B', -Q)
    return P
end

P = solve_riccati(A, B, Q, R)
println("Solution to the Riccati equation: ", P)

# Check controllability
function is_controllable(A, B)
    n = size(A, 1)
    ctrb_matrix = hcat([A^i * B for i in 0:n-1]...)
    rank(ctrb_matrix) == n
end

controllable = is_controllable(A, B)
println("Is the system controllable? ", controllable)