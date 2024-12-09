using LinearAlgebra, DifferentialEquations, ForwardDiff, KrylovKit
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

λ0 = 3.0 + 3.0im

u_end = find_steady_state(; P=P, ω_tilde, N_A=N_A, κ=κ, ω_r=ω_r, longmax=longmax, transmax=transmax, eigen_index=1, λ0=λ0)


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
eigenvalues = eigvals(A)
is_stable = all(real(eigenvalues) .< 0)
A = sparse(A)


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
function is_controllable(A, B; n=nothing)
    if n isa Nothing
        n = size(A, 1)
    end
    ctrb_matrix = hcat([A^i * B for i in 0:n-1]...)
    rank(ctrb_matrix) == n
end


control_vec = ComplexF64[P, ω_tilde, 0.0]
B = zeros(Float64, 2 * length(u_end), 2 * 3)
control_matrix!(B, complex_to_real(control_vec), nothing, nothing)
B = sparse(B)
controllable = is_controllable(A, B; n=40)


Q = Array{Float64}(I(2 * length(u_end)))
R = Array{Float64}(I(6))

function solve_riccati(A, B, Q, R)
    P = lyap(A, B * inv(R) * B', -Q)
    return P
end

P = solve_riccati(A, B, Q, R)
println("Solution to the Riccati equation: ", P)



println("Is the system controllable? ", controllable)