using DifferentialEquations, Plots, LaTeXStrings, LinearAlgebra, KrylovKit, NonlinearSolve, ForwardDiff
include("multimomenta_lib.jl")

P = 3.0
ω_tilde = 6.5
N_A = 10^5
κ = 8.1
E_0 = 40.0
ω_r = 0.05
longmax = 10
transmax = 10

ω_tildes = LinRange(6.5, 0.0, 50)

u_end, is_stable = find_steady_state_and_stability(; P=P, ω_tilde=ω_tildes[begin], N_A=N_A, κ=κ, ω_r=ω_r, longmax=longmax, transmax=transmax, eigen_index=1, λ0=3.0 + 3.0im);
λ0 = u_end[1]
λs = zeros(ComplexF64, length(ω_tildes))
markers = []

for (i, ω) in enumerate(ω_tildes)
    u_end, is_stable = find_steady_state_and_stability(; P=P, ω_tilde=ω, N_A=N_A, κ=κ, ω_r=ω_r, longmax=longmax, transmax=transmax, eigen_index=1, λ0=λ0)
    λ0 = u_end[1]
    λs[i] = λ0
    marker_style = is_stable ? :circle : :xcross
    push!(markers, marker_style)
end
scatter(real.(λs), imag.(λs), marker=markers, label="")
plot!(real.(λs[markers.==:circle]), imag.(λs[markers.==:circle]), seriestype=:scatter, marker=:circle, label="stable")
plot!(real.(λs[markers.==:xcross]), imag.(λs[markers.==:xcross]), seriestype=:scatter, marker=:xcross, label="unstable")

# Example usage
# f(z) = z^2
# plot_complex_function(f; xlims=(-3.0, 3.0), ylims=(-3.0, 3.0), gamma=0.5)

u_end, is_stable = find_steady_state_and_stability(; P=1.0, ω_tilde=4.0, N_A=N_A, κ=κ, ω_r=ω_r, longmax=longmax, transmax=transmax, eigen_index=1, λ0=3.0 + 3.0im);
plot_gaussian_bumps(u_end; longmax=longmax, transmax=transmax, pixel_per_bump=40, gamma=0.5, title="Normal Dicke Fixed Point")

u_end, is_stable = find_steady_state_and_stability(; P=0.6, ω_tilde=1.0, N_A=N_A, κ=κ, ω_r=ω_r, longmax=longmax, transmax=transmax, eigen_index=1, λ0=3.0im);
plot_gaussian_bumps(u_end; longmax=longmax, transmax=transmax, pixel_per_bump=40, gamma=0.5, title="Fixed Point in the Chaotic Region")