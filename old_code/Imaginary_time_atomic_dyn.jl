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

vec_dim = 1 + (2 * longmax + 1) * (2 * transmax + 1)
u0 = (1 - im) / (2 * sqrt(2 * N_A)) * ones(ComplexF64, vec_dim)
u0[1] = 1.0
u0[2+to_1d_index(0, 0, transmax, longmax)] = 1.0
u0[2+to_1d_index(1, 0, transmax, longmax)] = 0.0
u0[2+to_1d_index(0, 1, transmax, longmax)] = 0.0
u0[2+to_1d_index(-1, 0, transmax, longmax)] = 0.0
u0[2+to_1d_index(0, -1, transmax, longmax)] = 0.0
u0norm = sum(abs.(u0[2:end]) .^ 2)
u0[2:end] = u0[2:end] / sqrt(u0norm)

function drift!(du, u, p, t)
    u[2:end] = u[2:end] / norm(u[2:end])
    dispative_dynamics!(du, u, p, t; P=P, ω_tilde=ω_tilde, N_A=N_A, κ=κ, E_0=E_0, ω_r=ω_r, longmax=longmax, transmax=transmax)
    du[2:end] -= u[2:end] * dot(conj.(u[2:end]), du[2:end]) / dot(conj.(u[2:end]), u[2:end])
end

trecord = 0.0:0.05:5000.0
tspan = (trecord[begin], trecord[end])

probODE = ODEProblem(drift!, u0, tspan)
sol = solve(probODE, Tsit5(); reltol=10^-4, abstol=10^-4, dt=10^(-3), maxiters=10^12, save_noise=false, save_everystep=false, saveat=trecord)

p = plot(sol.t, map(x -> real(x[1]), sol.u), xlabel="Time", ylabel="Cavity Field", label=L"Re[$\lambda$]")
plot!(p, sol.t, map(x -> imag(x[1]), sol.u), xlabel="Time", ylabel="Cavity Field", label=L"Im[$\lambda$]")
savefig(p, "ImagTimeDynamics.png")

duTest = 0.0 * similar(u0)
multimomenta_model_drift!(duTest, sol.u[end], nothing, nothing; P=P, ω_tilde=ω_tilde, N_A=N_A, κ=κ, E_0=E_0, ω_r=ω_r, longmax=longmax, transmax=transmax)

println("norm of dlambda $(norm(duTest[1]))")

datom = duTest[2:end] - sol.u[end][2:end] * dot(conj.(sol.u[end][2:end]), duTest[2:end]) / dot(conj.(sol.u[end][2:end]), sol.u[end][2:end])

println("norm of datom $(norm(datom))")

dispative_dynamics!(duTest, sol.u[end], nothing, nothing; P=P, ω_tilde=ω_tilde, N_A=N_A, κ=κ, E_0=E_0, ω_r=ω_r, longmax=longmax, transmax=transmax)

println("norm of dlambda $(norm(duTest[1]))")

datom = duTest[2:end] - sol.u[end][2:end] * dot(conj.(sol.u[end][2:end]), duTest[2:end]) / dot(conj.(sol.u[end][2:end]), sol.u[end][2:end])

println("norm of datom $(norm(datom))")