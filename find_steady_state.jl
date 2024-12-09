using DifferentialEquations, Plots, LaTeXStrings, LinearAlgebra, KrylovKit, NonlinearSolve, ForwardDiff
include("multimomenta_lib.jl")

P = 0.5
ω_tilde = 2.5
N_A = 10^5
κ = 8.1
E_0 = 40.0
ω_r = 0.05
longmax = 10
transmax = 10

function find_steady_state_and_stability(; P=0.5, ω_tilde=1.5, N_A=10^5, κ=8.1, E_0=40.0, ω_r=0.05, longmax=2, transmax=2, eigen_index=1, λ0=3.0 + 3.0im)
    u_end = find_steady_state(; P=P, ω_tilde, N_A=N_A, κ=κ, ω_r=ω_r, longmax=longmax, transmax=transmax, eigen_index=eigen_index, λ0=λ0)

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
    println(eigenvalues[real(eigenvalues).>=10^(-4)])
    is_stable = all(real(eigenvalues) .< 10^(-3))
    return u_end, is_stable
end

ω_tildes = LinRange(2.5, 0.0, 50)

u_end, is_stable = find_steady_state_and_stability(; P=P, ω_tilde=ω_tildes[begin], N_A=N_A, κ=κ, ω_r=ω_r, longmax=longmax, transmax=transmax, eigen_index=2, λ0=3.0 + 3.0im);
λ0 = u_end[1]
λs = zeros(ComplexF64, length(ω_tildes))
markers = []

for (i, ω) in enumerate(ω_tildes)
    u_end, is_stable = find_steady_state_and_stability(; P=P, ω_tilde=ω, N_A=N_A, κ=κ, ω_r=ω_r, longmax=longmax, transmax=transmax, eigen_index=2, λ0=λ0)
    λ0 = u_end[1]
    λs[i] = λ0
    marker_style = is_stable ? :circle : :xcross
    push!(markers, marker_style)
end
scatter(real.(λs), imag.(λs), marker=markers, label="")
plot!(real.(λs[markers.==:circle]), imag.(λs[markers.==:circle]), seriestype=:scatter, marker=:circle, label="stable")
plot!(real.(λs[markers.==:xcross]), imag.(λs[markers.==:xcross]), seriestype=:scatter, marker=:xcross, label="unstable")

function plot_complex_function(f; xlims=(-1.0, 1.0), ylims=(-1.0, 1.0), resolution=100, gamma=1, kwargs...)
    x = range(xlims[1], xlims[2], length=resolution)
    y = range(ylims[1], ylims[2], length=resolution)
    z = [f(complex(re, im)) for re in x, im in y]

    plot_complex_mesh(x, y, z; xlims=xlims, ylims=ylims, gamma=gamma, color=:hsv, kwargs...)
end

function plot_complex_mesh(x, y, z; gamma=0.5, kwargs...)
    angles = angle.(z)
    mags = (abs.(z)) .^ gamma
    colors = HSV.((angles) ./ (π / 180.0), 1.0, mags ./ maximum(mags))

    heatmap(x, y, colors; kwargs...)#, kwargs...)
end

# Example usage
# f(z) = z^2
# plot_complex_function(f; xlims=(-3.0, 3.0), ylims=(-3.0, 3.0), gamma=0.5)

function plot_gaussian_bumps(u0; longmax=2, transmax=2, sigma=0.1, pixel_per_bump=10, kwargs...)
    x = range(-(transmax + 1), (transmax + 1), length=pixel_per_bump * (2 * transmax + 1))
    y = range(-(longmax + 1), (longmax + 1), length=pixel_per_bump * (2 * longmax + 1))
    X, Y = [x for x in x, y in y], [y for x in x, y in y]
    Z = zeros(size(X))

    Threads.@threads for n in -transmax:transmax
        for m in -longmax:longmax
            idx = to_1d_index(n, m, transmax, longmax)
            amplitude = abs(u0[2+idx])
            Z .+= amplitude * exp.(-((X .- n) .^ 2 .+ (Y .- m) .^ 2) / (2 * sigma^2))
        end
    end

    plot_complex_mesh(x, y, Z; color=:viridis, xlabel=L"$\frac{k_x}{k_r}$", ylabel=L"$\frac{k_z}{k_r}$", xlims=(x[begin], x[end]), ylim=(y[begin], y[end]), kwargs...)
end

# Example usage
plot_gaussian_bumps(u_end; longmax=longmax, transmax=transmax, pixel_per_bump=40, gamma=0.5)

