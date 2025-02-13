using LinearAlgebra, DifferentialEquations, ForwardDiff, KrylovKit, Plots, LaTeXStrings
include("multimomenta_lib.jl")

P = 1.0
ω_tilde = 6.5
N_A = 10^5
κ = 100.0
E_0 = 40.0
ω_r = 0.05
longmax = 2
transmax = 2
ω_tilde = 6.5

λ0 = 0.562 + 0.0215im

u_end = find_steady_state(; P=P, ω_tilde, N_A=N_A, κ=κ, ω_r=ω_r, longmax=longmax, transmax=transmax, eigen_index=1, λ0=λ0)
du = similar(u_end)
du .= 0.0
multimomenta_model_drift!(du, u_end, nothing, nothing; P=P, ω_tilde, N_A=N_A, κ=κ, ω_r=ω_r, longmax=longmax, transmax=transmax)
@assert norm(du[1]) < 10^-6
@assert norm(du[2:end] - (u_end[2:end]' * du[2:end]) / (u_end[2:end]' * u_end[2:end]) * u_end[2:end]) < 10^-9

vec_dim = 1 + (2 * longmax + 1) * (2 * transmax + 1)
u0 = (1 - im) / (2 * sqrt(2 * N_A)) * ones(ComplexF64, vec_dim)
u0[1] = -0.01
u0[2+to_1d_index(0, 0, transmax, longmax)] = 1.0
u0[2+to_1d_index(1, 0, transmax, longmax)] = 0.0
u0[2+to_1d_index(0, 1, transmax, longmax)] = 0.0
u0[2+to_1d_index(-1, 0, transmax, longmax)] = 0.0
u0[2+to_1d_index(0, -1, transmax, longmax)] = 0.0
u0[2:end] /= norm(u0[2:end])

function smoothstep(t)
    if t < 0
        return 0.0
    elseif t < 1
        return 3 * t^2 - 2 * t^3
    else
        return 1.0
    end
end

function jacobian_multimomenta_model_drift!(J, u, p, t)
    J .= ForwardDiff.jacobian(temp_u -> begin
            uComplex = real_to_complex(temp_u)
            du = similar(uComplex)
            multimomenta_model_drift!(du, uComplex, p, t; P=P, ω_tilde=ω_tilde, N_A=N_A, κ=κ, E_0=E_0, ω_r=ω_r, longmax=longmax, transmax=transmax)
            du[2:end] -= (uComplex[2:end]' * du[2:end]) / (uComplex[2:end]' * uComplex[2:end]) * uComplex[2:end]
            return complex_to_real(du)
        end, u)
end

function control_matrix!(B, st_control, p, t)
    B .= ForwardDiff.jacobian(control -> begin
            du = similar(u_end, typeof(complex(control[1], control[2])))
            u_end[2:end] /= norm(u_end[2:end])
            multimomenta_model_drift!(du, u_end, p, t; P=complex(control[1], control[2]), ω_tilde=control[3], N_A=N_A, κ=κ, E_0=E_0, ω_r=ω_r, ϵ=complex(control[4], control[5]), longmax=longmax, transmax=transmax)
            du[2:end] -= (u_end[2:end]' * du[2:end]) / (u_end[2:end]' * u_end[2:end]) * u_end[2:end]
            return complex_to_real(du)
        end, st_control)
end

function solve_riccati(A, B, Q, R)
    AREH = [A -B*inv(R)*B'; -Q -A']
    evals, evecs = eigen(AREH)
    n = size(A, 1)
    return evecs[n+1:2n, 1:n] * inv(evecs[1:n, 1:n])
end

function dichotomous_subspace_projector(A::Matrix{T}; tol=1e-6)::Matrix{T} where {T<:Real}
    evals, evecs = eigen(A)
    non_marg_osc_indices = findall(x -> (abs(real(x)) > tol) && imag(x) >= 0, evals)
    proj = hcat(real.(evecs[:, non_marg_osc_indices]), imag.(evecs[:, non_marg_osc_indices]))
    proj = remove_zero_columns!(proj)
    for i in size(proj, 2)
        proj[:, i] /= norm(proj[:, i])
    end
    U, s, Vt = svd(proj)
    return U * Vt
end

function sde_drift!(du, u, p, t)
    # Update A and B matrices
    A = zeros(Float64, 2 * length(u), 2 * length(u))
    jacobian_multimomenta_model_drift!(A, complex_to_real(u), p, t)
    proj = dichotomous_subspace_projector(A; tol=1e-8)
    Adc = proj' * A * proj

    control_vec = Float64[P, 0.0, ω_tilde, 0.0, 0.0]
    B = zeros(Float64, 2 * length(u), length(control_vec))
    control_matrix!(B, control_vec, p, t)
    Bdc = proj' * B

    Q = zeros(Float64, size(A))
    Q[1, 1] = 1.0
    Qdc = proj' * Q * proj
    R = Array{Float64}(I(length(control_vec))) * 100

    Ricca = solve_riccati(Array(Adc), Array(Bdc), Qdc, R)
    CostMat = real.(Ricca)
    L = inv(R) * Bdc' * CostMat * proj'

    cont_vec = -L * complex_to_real(u - u_end)
    multimomenta_model_drift!(du, u, p, t; P=P * smoothstep(t / 600.0) + complex(cont_vec[1], cont_vec[2]), ω_tilde=ω_tilde + cont_vec[3], N_A=N_A, κ=κ, E_0=E_0, ω_r=ω_r, ϵ=complex(cont_vec[4], cont_vec[5]), longmax=longmax, transmax=transmax)
end

trecord = 0.0:1.0:100.0
tspan = (trecord[begin], trecord[end])

u0 = find_steady_state(; P=P, ω_tilde, N_A=N_A, κ=κ, ω_r=ω_r, longmax=longmax, transmax=transmax, eigen_index=1, λ0=λ0)
du = similar(u_end)
du .= 0.0
multimomenta_model_drift!(du, u_end, nothing, nothing; P=P, ω_tilde, N_A=N_A, κ=κ, ω_r=ω_r, longmax=longmax, transmax=transmax)
@assert norm(du[1]) < 10^-6
@assert norm(du[2:end] - (u_end[2:end]' * du[2:end]) / (u_end[2:end]' * u_end[2:end]) * u_end[2:end]) < 10^-9
u0[1] += 0.001

probODE = ODEProblem(sde_drift!, u0, tspan)
sol = solve(probODE, Tsit5(); reltol=10^-4, abstol=10^-4, dt=10^(-3), maxiters=10^7, save_noise=false, save_everystep=false, saveat=trecord)

p = plot(sol.t, map(x -> real(x[1]), sol.u))
savefig(p, "state_dependent_lqg.png")

p = scatter(map(x -> real(x[1]), sol.u), map(x -> imag(x[1]), sol.u), zcolor=sol.t, xlims=(-8, 8), ylims=(-8, 8), colormap=:viridis, colorbar=true, markerstrokewidth=0, markersize=0.8, label="P=$(P), ω=$(ω_tilde)", dpi=300)
scatter!(p, [-real(u_end[1])], [-imag(u_end[1])], label="Target State")
title!(p, "Trajectory with State-Dependent LQG Control")
xaxis!(p, L"Re[\lambda]")
yaxis!(p, L"Im[\lambda]")
savefig(p, "state_dependent_lqg_trajectory.png")