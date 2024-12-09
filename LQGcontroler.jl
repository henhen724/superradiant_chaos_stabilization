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


control_vec = ComplexF64[P, ω_tilde, 0.0]
B = zeros(Float64, 2 * length(u_end), 2 * length(control_vec))
control_matrix!(B, complex_to_real(control_vec), nothing, nothing)
B = sparse(B)
# Proj = controllable_subspace_projector(A, B; tol=1e-10)
# controllable = controllable(A, B)

Q = Array{Float64}(I(2 * length(u_end)))
R = Array{Float64}(I(2 * length(control_vec))) / 3

function solve_riccati(A, B, Q, R)
    AREH = [A -B*inv(R)*B'; -Q -A']
    evals, evecs = eigen(AREH)
    n = size(A, 1)
    println(evals[1:n])
    return evecs[n+1:2n, 1:n] * inv(evecs[1:n, 1:n])
end

CostMat = solve_riccati(Array(A), Array(B), Q, R)


vec_dim = 1 + (2 * longmax + 1) * (2 * transmax + 1)
u0 = (1 - im) / (2 * sqrt(2 * N_A)) * ones(ComplexF64, vec_dim)
u0[2+to_1d_index(0, 0, transmax, longmax)] = 1.0
u0[2+to_1d_index(1, 0, transmax, longmax)] = 0.0
u0[2+to_1d_index(0, 1, transmax, longmax)] = 0.0
u0[2+to_1d_index(-1, 0, transmax, longmax)] = 0.0
u0[2+to_1d_index(0, -1, transmax, longmax)] = 0.0
u0norm = sum(abs.(u0[2:end]) .^ 2)
u0[2:end] = u0[2:end] / sqrt(u0norm)

function sde_drift!(du, u, p, t)
    cont_vec = real_to_complex(inv(R) * B' * CostMat * complex_to_real(u - u_end))
    multimomenta_model_drift!(du, u, p, t; P=P + cont_vec[1], ω_tilde=ω_tilde + cont_vec[2], N_A=N_A, κ=κ, E_0=E_0, ω_r=ω_r, ϵ=cont_vec[3], longmax=longmax, transmax=transmax)
end

trecord = 0.0:0.05:5000.0
tspan = (trecord[begin], trecord[end])

probODE = ODEProblem(sde_drift!, u0, tspan)
# sol = solve(probSDE, RKMilGeneral(; ii_approx=IICommutative()); adaptive=false, dt=2^(-15))
sol = solve(probODE, Tsit5(); reltol=10^-4, abstol=10^-4, dt=10^(-3), maxiters=10^12, save_noise=false, save_everystep=false, saveat=trecord)
