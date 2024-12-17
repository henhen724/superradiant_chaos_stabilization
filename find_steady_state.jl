using DifferentialEquations, Plots, LaTeXStrings, LinearAlgebra, KrylovKit, NonlinearSolve, ForwardDiff
include("multimomenta_lib.jl")

P = 1.0
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
largest_eigs = zeros(ComplexF64, length(ω_tildes))
markers = []
u_ends = []

for (i, ω) in enumerate(ω_tildes)
    u_end, is_stable, largest_eigs[i] = find_steady_state_and_stability(; P=P, ω_tilde=ω, N_A=N_A, κ=κ, ω_r=ω_r, longmax=longmax, transmax=transmax, eigen_index=1, λ0=λ0)
    du = similar(u_end)
    du .= 0.0
    multimomenta_model_drift!(du, u_end, nothing, nothing; P=P, ω_tilde=ω, N_A=N_A, κ=κ, ω_r=ω_r, longmax=longmax, transmax=transmax)
    println(du[1])
    @assert norm(du[1]) < 10^-6
    if abs(u_end[2+to_1d_index(transmax, longmax, transmax, longmax)]) > 10^-3
        println("Warning momentum cutoff may effect the result")
        plot_gaussian_bumps(u_end; longmax=longmax, transmax=transmax, pixel_per_bump=40)
        sleep(0.5)
    end
    # println(norm(du[2:end] - (u_end[2:end]' * du[2:end]) / (u_end[2:end]' * u_end[2:end]) * u_end[2:end]))
    @assert norm(du[2:end] - (u_end[2:end]' * du[2:end]) / (u_end[2:end]' * u_end[2:end]) * u_end[2:end]) < 10^-6
    λ0 = u_end[1]
    λs[i] = λ0
    marker_style = is_stable ? :circle : :xcross
    push!(markers, marker_style)
    push!(u_ends, u_end)
end
p = scatter(real.(λs), imag.(λs), zcolor=ω_tildes, marker=markers, label="", dpi=400)
plot!(p, real.(λs[markers.==:circle]), imag.(λs[markers.==:circle]), seriestype=:scatter, marker=:circle, label="stable")
plot!(p, real.(λs[markers.==:xcross]), imag.(λs[markers.==:xcross]), seriestype=:scatter, marker=:xcross, label="unstable")
xlabel!(p, L"Re[λ]")
ylabel!(p, L"Im[λ]")
title!(p, "Cavity Field at Fixed Points")
savefig(p, "plot.png")

index = 48
p = plot_gaussian_bumps(u_ends[index]; transmax=transmax, longmax=longmax, dpi=400, pixel_per_bump=40)
annotate!(p, -18.0, -10.0, text("\$\\tilde{\\omega}\$=$(round(ω_tildes[index],digits=1))", :left, 12; color=:black))
savefig(p, "plot.png")

p = scatter(repeat([P], length(ω_tildes)), ω_tildes, xlims=(0.0, 1.25), ylims=(-0.1, 7.0), label="", marker=markers, left_margin=3Plots.mm, dpi=400)
scatter(repeat([P], length(ω_tildes)), ω_tildes)
plot!(p, repeat([P], length(ω_tildes)), ω_tildes[markers.==:circle], seriestype=:scatter, marker=:circle, label="stable")
plot!(p, repeat([P], length(ω_tildes)), ω_tildes[markers.==:xcross], seriestype=:scatter, marker=:xcross, label="unstable")
title!(p, "Phase Diagram")
xlabel!(p, L"P=\frac{\Omega}{\sqrt{4 \Delta_a \omega_r}}")
ylabel!(p, L"\tilde{ω}=\frac{\omega_c}{E_0} = \frac{4 \Delta_a \omega_c}{g^2_0 N}")
hline!(p, [2.0], label="Blue Dispersive Shift", color=:blue)
savefig(p, "plot.png")

# Example usage
# f(z) = z^2
# plot_complex_function(f; xlims=(-3.0, 3.0), ylims=(-3.0, 3.0), gamma=0.5)

u_end, is_stable = find_steady_state_and_stability(; P=1.0, ω_tilde=4.0, N_A=N_A, κ=κ, ω_r=ω_r, longmax=longmax, transmax=transmax, eigen_index=1, λ0=3.0 + 3.0im);
plot_gaussian_bumps(u_end; longmax=longmax, transmax=transmax, pixel_per_bump=40, gamma=0.5, title="Normal Dicke Fixed Point")

u_end, is_stable = find_steady_state_and_stability(; P=0.6, ω_tilde=1.0, N_A=N_A, κ=κ, ω_r=ω_r, longmax=longmax, transmax=transmax, eigen_index=1, λ0=3.0im);
plot_gaussian_bumps(u_end; longmax=longmax, transmax=transmax, pixel_per_bump=40, gamma=0.5, title="Fixed Point in the Chaotic Region")