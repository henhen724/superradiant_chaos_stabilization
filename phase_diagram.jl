using DifferentialEquations, Plots, LaTeXStrings, LinearAlgebra, KrylovKit, NonlinearSolve, ForwardDiff
include("multimomenta_lib.jl")

N_A = 10^5
κ = 8.1
E_0 = 40.0
ω_r = 0.05
longmax = 10
transmax = 10

Ps = LinRange(0.0, 1.0, 13)
ω_tildes = LinRange(0.0, 6.0, 13)

function scan_lambda0(P, ω_tilde; N_A=10^5, κ=8.1, ω_r=0.05, longmax=10, transmax=10)
    λ0_real_values = LinRange(-3.0, 3.0, 5)
    λ0_imag_values = LinRange(-3.0, 3.0, 5)
    λ0_values = [re_p + im_p * im for re_p in λ0_real_values, im_p in λ0_imag_values]
    results = []
    for λ0 in λ0_values
        u_end, is_stable = find_steady_state_and_stability(; P=P, ω_tilde=ω_tilde, N_A=N_A, κ=κ, ω_r=ω_r, longmax=longmax, transmax=transmax, eigen_index=1, λ0=λ0)
        if all(abs(u_end[1] - res[1]) .> 0.01 for res in results)
            push!(results, (u_end[1], is_stable))
        end
    end
    return results
end

@enum Phase Normal Stable_SR Chaotic Bistable Error

function classify_phase(results)::Phase
    zero_stable = any(abs(res[1]) < 0.01 && res[2] for res in results)
    other_stable = any(abs(res[1]) > 0.01 && res[2] for res in results)

    if zero_stable
        if other_stable
            return Bistable
        else
            return Normal
        end
    else
        if other_stable
            return Stable_SR
        else
            return Chaotic
        end
    end
end

phases = Matrix{Union{Nothing,Phase}}(nothing, length(Ps), length(ω_tildes))

Threads.@threads for i in 1:length(Ps)
    P = Ps[i]
    for j in 1:length(ω_tildes)
        ω_tilde = ω_tildes[j]
        try
            results = scan_lambda0(P, ω_tilde)
            phases[i, j] = classify_phase(results)
        catch e
            println("Error encountered: ", e)
            results = []
            phases[i, j] = classify_phase(results)
        end
    end
end



p = heatmap(Ps, ω_tildes, phases, c=:viridis, xlabel="P", ylabel=L"\tilde{\omega}", title="Phase Diagram", colorbar=false)
