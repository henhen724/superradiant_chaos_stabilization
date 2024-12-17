using LinearAlgebra, DifferentialEquations, ForwardDiff, KrylovKit, Plots, LaTeXStrings
include("multimomenta_lib.jl")

P = 1.0
ω_tilde = 6.5
N_A = 10^5
κ = 8.1
E_0 = 40.0
ω_r = 0.05
longmax = 10
transmax = 10
ω_tilde = 6.5

λ0 = 0.562 + 0.0215im

u_end = find_steady_state(; P=P, ω_tilde, N_A=N_A, κ=κ, ω_r=ω_r, longmax=longmax, transmax=transmax, eigen_index=1, λ0=λ0)
du = similar(u_end)
du .= 0.0
multimomenta_model_drift!(du, u_end, nothing, nothing; P=P, ω_tilde, N_A=N_A, κ=κ, ω_r=ω_r, longmax=longmax, transmax=transmax);
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

function sde_drift!(du, u, p, t)
    # cont_vec = -L * complex_to_real(u - u_end)
    cont_vec = zeros(5)
    # println(cont_vec)
    multimomenta_model_drift!(du, u, p, t; P=P * smoothstep(t / 600.0) + complex(cont_vec[1], cont_vec[2]), ω_tilde=ω_tilde + cont_vec[3], N_A=N_A, κ=κ, E_0=E_0, ω_r=ω_r, ϵ=complex(cont_vec[4], cont_vec[5]), longmax=longmax, transmax=transmax)
end

trecord = 0.0:0.5:5000.0
tspan = (trecord[begin], trecord[end])

probODE = ODEProblem(sde_drift!, u0, tspan)
# sol = solve(probSDE, RKMilGeneral(; ii_approx=IICommutative()); adaptive=false, dt=2^(-15))
sol = solve(probODE, Tsit5(); reltol=10^-4, abstol=10^-4, dt=10^(-3), maxiters=10^12, save_noise=false, save_everystep=false, saveat=trecord)

p = plot(sol.t, map(x -> real(x[1]), sol.u))
savefig(p, "weird_meta_stable.png")

p = scatter(map(x -> real(x[1]), sol.u), map(x -> imag(x[1]), sol.u), zcolor=sol.t, xlims=(-8, 8), ylims=(-8, 8), colormap=:viridis, colorbar=true, markerstrokewidth=0, markersize=0.8, label="P=$(P), ω=$(ω_tilde)", dpi=300)
scatter!(p, [-real(u_end[1])], [-imag(u_end[1])], label="Target State")
title!(p, "Trajectory without Heterodyne Measurement")
xaxis!(p, L"Re[\lambda]")
yaxis!(p, L"Im[\lambda]")
savefig(p, "weird_meta_stable.png")

u_normal, is_normal_stable = find_steady_state_and_stability(; P=P, ω_tilde, N_A=N_A, κ=κ, ω_r=ω_r, longmax=longmax, transmax=transmax, eigen_index=1, λ0=0.0 + 0.0im)
println("Is normal stable? ", is_normal_stable)

# Define Jacobian function
function jacobian_multimomenta_model_drift!(J, u, p, t)
    J .= ForwardDiff.jacobian(temp_u -> begin
            uComplex = real_to_complex(temp_u)
            du = similar(uComplex)
            multimomenta_model_drift!(du, uComplex, p, t; P=P, ω_tilde=ω_tilde, N_A=N_A, κ=κ, E_0=E_0, ω_r=ω_r, longmax=longmax, transmax=transmax)
            du[2:end] -= (uComplex[2:end]' * du[2:end]) / (uComplex[2:end]' * uComplex[2:end]) * uComplex[2:end]
            return complex_to_real(du)
        end, u)
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

# Calculate Jacobian at u_end
A = zeros(Float64, 2 * length(u_end), 2 * length(u_end))
jacobian_multimomenta_model_drift!(A, complex_to_real(u_end), nothing, 0.0)
evals, evecs = eigen(A)
is_stable = all(real(evals) .< 10^-6)
proj = dichotomous_subspace_projector(A; tol=1e-8)

Adc = proj' * A * proj

println("Is the fixed point stable? ", is_stable)

# Define control matrix function
function control_matrix!(B, st_control, p, t)
    B .= ForwardDiff.jacobian(control -> begin
            # cont_complex = real_to_complex(control)
            du = similar(u_end, typeof(complex(control[1], control[2])))
            u_end[2:end] /= norm(u_end[2:end])
            multimomenta_model_drift!(du, u_end, p, t; P=complex(control[1], control[2]), ω_tilde=control[3], N_A=N_A, κ=κ, E_0=E_0, ω_r=ω_r, ϵ=complex(control[4], control[5]), longmax=longmax, transmax=transmax)
            du[2:end] -= (u_end[2:end]' * du[2:end]) / (u_end[2:end]' * u_end[2:end]) * u_end[2:end]
            return complex_to_real(du)
        end, st_control)
end

control_vec = Float64[P, 0.0, ω_tilde, 0.0, 0.0]
B = zeros(Float64, 2 * length(u_end), length(control_vec))
control_matrix!(B, control_vec, nothing, nothing)
Bdc = proj' * B
# Proj = controllable_subspace_projector(A, B; tol=1e-10)
# controllable = controllable(A, B)

Q = zeros(Float64, size(A))
Q[1, 1] = 1.0
Qdc = proj' * Q * proj
R = Array{Float64}(I(length(control_vec))) * 100


function solve_riccati(A, B, Q, R)
    AREH = [A -B*inv(R)*B'; -Q -A']
    evals, evecs = eigen(AREH)
    n = size(A, 1)
    println(evals[1:n])
    return evecs[n+1:2n, 1:n] * inv(evecs[1:n, 1:n])
end

Ricca = solve_riccati(Array(Adc), Array(Bdc), Qdc, R)
println(norm(imag.(Ricca)))
CostMat = real.(Ricca)
println(norm(real.(Ricca)))
L = inv(R) * Bdc' * CostMat * proj'
println(norm(L))

function smoothstep(t)
    if t < 0
        return 0.0
    elseif t < 1
        return 3 * t^2 - 2 * t^3
    else
        return 1.0
    end
end

function sde_drift!(du, u, p, t)
    # if t > 1000.0
    #     cont_vec = -L * complex_to_real(u - u_end)
    # else
    #     cont_vec = zeros(5)
    # end
    cont_vec = -L * complex_to_real(u - u_end)
    # println(cont_vec)
    multimomenta_model_drift!(du, u, p, t; P=P * smoothstep(t / 600.0) + complex(cont_vec[1], cont_vec[2]), ω_tilde=ω_tilde + cont_vec[3], N_A=N_A, κ=κ, E_0=E_0, ω_r=ω_r, ϵ=complex(cont_vec[4], cont_vec[5]), longmax=longmax, transmax=transmax)
end

trecord = 0.0:0.5:1500.0
tspan = (trecord[begin], trecord[end])

u0 = find_steady_state(; P=P, ω_tilde, N_A=N_A, κ=κ, ω_r=ω_r, longmax=longmax, transmax=transmax, eigen_index=1, λ0=λ0)
du = similar(u_end)
du .= 0.0
multimomenta_model_drift!(du, u_end, nothing, nothing; P=P, ω_tilde, N_A=N_A, κ=κ, ω_r=ω_r, longmax=longmax, transmax=transmax);
@assert norm(du[1]) < 10^-6
@assert norm(du[2:end] - (u_end[2:end]' * du[2:end]) / (u_end[2:end]' * u_end[2:end]) * u_end[2:end]) < 10^-9
u0[1] += 0.001

probODE = ODEProblem(sde_drift!, u0, tspan)
# sol = solve(probSDE, RKMilGeneral(; ii_approx=IICommutative()); adaptive=false, dt=2^(-15))
sol = solve(probODE, Tsit5(); reltol=10^-4, abstol=10^-4, dt=10^(-3), maxiters=10^7, save_noise=false, save_everystep=false, saveat=trecord)

p = plot(sol.t, map(x -> real(x[1]), sol.u))

cont_vec = zeros(5)

function open_loop_sde_drift!(du, u, p, t)
    # if t > 1000.0
    #     cont_vec = -L * complex_to_real(u - u_end)
    # else
    #     cont_vec = zeros(5)
    # end
    # println(cont_vec)
    multimomenta_model_drift!(du, u, p, t; P=P * smoothstep(t / 600.0) + complex(cont_vec[1], cont_vec[2]), ω_tilde=ω_tilde + cont_vec[3], N_A=N_A, κ=κ, E_0=E_0, ω_r=ω_r, ϵ=complex(cont_vec[4], cont_vec[5]), longmax=longmax, transmax=transmax)
end

probODE = ODEProblem(sde_drift!, u0, tspan)
# sol = solve(probSDE, RKMilGeneral(; ii_approx=IICommutative()); adaptive=false, dt=2^(-15))
sol_open = solve(probODE, Tsit5(); reltol=10^-4, abstol=10^-4, dt=10^(-3), maxiters=10^7, save_noise=false, save_everystep=false, saveat=trecord)

p = plot(sol.t, map(x -> real(x[1]), sol_open.u))

plot_gaussian_bumps(u_end; longmax=longmax, transmax=transmax, pixel_per_bump=40, gamma=0.5, title="Normal Dicke Fixed Point")
plot_gaussian_bumps(sol.u[end]; longmax=longmax, transmax=transmax, pixel_per_bump=40, gamma=0.5, title="Normal Dicke Fixed Point")