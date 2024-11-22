using DifferentialEquations, Plots, LaTeXStrings

function to_1d_index(n::Int, m::Int, transmax::Int, longmax::Int)::Int
    @assert -transmax <= n <= transmax
    @assert -longmax <= m <= longmax
    return (n + transmax) * (2 * longmax + 1) + m + longmax
end

function safe_index_2D(vec::Vector{T}, n::Int, m::Int, transmax::Int, longmax::Int) where {T}
    if -transmax <= n <= transmax && -longmax <= m <= longmax
        return vec[2+to_1d_index(n, m, transmax, longmax)]
    else
        return zero(T)
    end
end

function run_dynamics_with_hetero(P, ω_tilde; N_A=10^5, κ=8.1, E_0=40.0, ω_r=0.05, longmax=10, transmax=10)
    ω_c = ω_tilde * E_0


    vec_dim = 1 + (2 * longmax + 1) * (2 * transmax + 1)
    u0 = (1 - im) / (2 * sqrt(2 * N_A)) * ones(ComplexF64, vec_dim)
    u0[2+to_1d_index(0, 0, transmax, longmax)] = 1.0
    u0[2+to_1d_index(1, 0, transmax, longmax)] = 0.0
    u0[2+to_1d_index(0, 1, transmax, longmax)] = 0.0
    u0[2+to_1d_index(-1, 0, transmax, longmax)] = 0.0
    u0[2+to_1d_index(0, -1, transmax, longmax)] = 0.0
    norm = sum(abs.(u0[2:end]) .^ 2)
    u0[2:end] = u0[2:end] / sqrt(norm)

    noise_prototype = zeros(ComplexF64, (vec_dim, 2)) #hetrodyne
    noise = RealWienerProcess!(0.0, zeros(2), zeros(2), save_everystep=false)
    # noise_prototype = zeros(Float64, (vec_dim, 1)) #homodyne

    function sde_drift(du, u, p, t)
        long_sum_all = 0
        checker_board_all = 0
        for n in -transmax:transmax
            for m in -longmax:longmax
                mom_indx = to_1d_index(n, m, transmax, longmax)
                trans_sum = safe_index_2D(u, n + 2, m, longmax, transmax) + safe_index_2D(u, n - 2, m, longmax, transmax)
                long_sum = safe_index_2D(u, n, m + 2, longmax, transmax) + safe_index_2D(u, n, m - 2, longmax, transmax)
                long_sum_all += conj(long_sum) * u[2+mom_indx]
                checker_board = safe_index_2D(u, n + 1, m + 1, longmax, transmax) +
                                safe_index_2D(u, n + 1, m - 1, longmax, transmax) +
                                safe_index_2D(u, n - 1, m + 1, longmax, transmax) +
                                safe_index_2D(u, n - 1, m - 1, longmax, transmax)
                checker_board_all += conj(checker_board) * u[2+mom_indx]
                du[2+mom_indx] = -im * ω_r * ((n^2 + m^2) * u[2+mom_indx] - conj(u[1]) * u[1] * long_sum - P * (u[1] + conj(u[1])) * checker_board - P^2 * trans_sum)
            end
        end
        du[1] = -(κ + im * ω_c) * u[1] + im * E_0 * u[1] * long_sum_all + im * E_0 * P * checker_board_all
    end

    function sde_diffusion(du, u, p, t)
        du[1, 1] = sqrt(E_0 * κ / (4 * N_A * ω_r))
        du[1, 2] = im * sqrt(E_0 * κ / (4 * N_A * ω_r)) # comment out for homodyne
    end

    trecord = 0.0:0.05:5000.0
    tspan = (trecord[begin], trecord[end])

    probSDE = SDEProblem(sde_drift, sde_diffusion, u0, tspan, noise_rate_prototype=noise_prototype, noise=noise, save_noise=false)
    # probODE = ODEProblem(sde_drift, u0, tspan)
    # sol = solve(probSDE, RKMilGeneral(; ii_approx=IICommutative()); adaptive=false, dt=2^(-15))
    sol = solve(probSDE, SOSRA2(); reltol=10^-3, abstol=10^-3, dt=10^(-3), maxiters=10^13, save_noise=false, save_everystep=false, saveat=trecord, dtmin=10^-15)
    return sol
end

function run_dynamics_no_meas(P, ω_tilde; N_A=10^5, κ=8.1, E_0=40.0, ω_r=0.05, longmax=10, transmax=10)
    ω_c = ω_tilde * E_0


    vec_dim = 1 + (2 * longmax + 1) * (2 * transmax + 1)
    u0 = (1 - im) / (2 * sqrt(2 * N_A)) * ones(ComplexF64, vec_dim)
    u0[2+to_1d_index(0, 0, transmax, longmax)] = 1.0
    u0[2+to_1d_index(1, 0, transmax, longmax)] = 0.0
    u0[2+to_1d_index(0, 1, transmax, longmax)] = 0.0
    u0[2+to_1d_index(-1, 0, transmax, longmax)] = 0.0
    u0[2+to_1d_index(0, -1, transmax, longmax)] = 0.0
    norm = sum(abs.(u0[2:end]) .^ 2)
    u0[2:end] = u0[2:end] / sqrt(norm)

    function sde_drift(du, u, p, t)
        long_sum_all = 0
        checker_board_all = 0
        for n in -transmax:transmax
            for m in -longmax:longmax
                mom_indx = to_1d_index(n, m, transmax, longmax)
                trans_sum = safe_index_2D(u, n + 2, m, longmax, transmax) + safe_index_2D(u, n - 2, m, longmax, transmax)
                long_sum = safe_index_2D(u, n, m + 2, longmax, transmax) + safe_index_2D(u, n, m - 2, longmax, transmax)
                long_sum_all += conj(long_sum) * u[2+mom_indx]
                checker_board = safe_index_2D(u, n + 1, m + 1, longmax, transmax) +
                                safe_index_2D(u, n + 1, m - 1, longmax, transmax) +
                                safe_index_2D(u, n - 1, m + 1, longmax, transmax) +
                                safe_index_2D(u, n - 1, m - 1, longmax, transmax)
                checker_board_all += conj(checker_board) * u[2+mom_indx]
                du[2+mom_indx] = -im * ω_r * ((n^2 + m^2) * u[2+mom_indx] - conj(u[1]) * u[1] * long_sum - P * (u[1] + conj(u[1])) * checker_board - P^2 * trans_sum)
            end
        end
        du[1] = -(κ + im * ω_c) * u[1] + im * E_0 * u[1] * long_sum_all + im * E_0 * P * checker_board_all
    end

    trecord = 0.0:0.05:5000.0
    tspan = (trecord[begin], trecord[end])

    probODE = ODEProblem(sde_drift, u0, tspan)
    # sol = solve(probSDE, RKMilGeneral(; ii_approx=IICommutative()); adaptive=false, dt=2^(-15))
    sol = solve(probODE, Tsit5(); reltol=10^-4, abstol=10^-4, dt=10^(-3), maxiters=10^12, save_noise=false, save_everystep=false, saveat=trecord)
    return sol
end

Ps = [0.238, 0.367, 0.450, 0.533, 0.617, 0.700]
ω_tildes = [1.880, 1.520, 1.170, 0.825, 0.475, 0.125]
P = Ps[4]
ω_tilde = ω_tildes[2]
sol = run_dynamics_with_hetero(P, ω_tilde; N_A=10^5)
plot(sol.t, map(x -> real(x[1]), sol.u))
p33 = plot(sol.t, map(x -> real(x[2+to_1d_index(3, 3, 10, 10)]), sol.u), label=L"Re[\langle c_{3,3} \rangle]")
plot!(p33, sol.t, map(x -> imag(x[2+to_1d_index(3, 3, 10, 10)]), sol.u), label=L"Im[\langle c_{3,3} \rangle]")
xaxis!(p33, L"Time [\mu s]")
yaxis!(p33, (-0.2, 0.2))
title!(p33, L"3,3 $k_r$ state $N_A=10^5$")
savefig(p33, "state33_timeseries.png")


p = scatter(map(x -> real(x[1]), sol.u), map(x -> imag(x[1]), sol.u), zcolor=sol.t, xlims=(-8, 8), ylims=(-8, 8), colormap=:viridis, colorbar=true, markerstrokewidth=0, markersize=0.8, label="P=$(P), ω=$(ω_tilde)")
title!(p, "Trajectory with Heterodyne Measurement")
xaxis!(L"Re[\lambda]")
yaxis!(L"Im[\lambda]")
savefig(p, "buildup_hetero.png")

sol_no_meas = run_dynamics_no_meas(P, ω_tilde)
plot(sol.t, map(x -> real(x[2+to_1d_index(3, 3, 10, 10)]), sol_no_meas.u))
ptimeseries = plot(sol.t, map(x -> real(x[1]), sol.u), label=L"Re[$\langle a \rangle$]")
plot!(ptimeseries, sol.t, map(x -> imag(x[1]), sol.u), label=L"Im[$\langle a \rangle$]")
title!(ptimeseries, "Time Series of Meanfield Dyanmics")
xaxis!(ptimeseries, L"Time [$\mu$s]")
savefig(ptimeseries, "traj_timeseries.png")
p = scatter(map(x -> real(x[1]), sol.u), map(x -> imag(x[1]), sol.u), zcolor=sol.t, xlims=(-8, 8), ylims=(-8, 8), colormap=:viridis, colorbar=true, markerstrokewidth=0, markersize=0.8, label="P=$(P), ω=$(ω_tilde)")
title!(p, "Trajectory without Measurement")
xaxis!(L"Re[\lambda]")
yaxis!(L"Im[\lambda]")
savefig(p, "traj_buildup.png")

N_As = [10^6, 10^5, 10^4, 10^3, 10^2, 10^1]
plot_list = []
for N_A in N_As
    sol = run_dynamics_with_hetero(P, ω_tilde; N_A=N_A)
    p = scatter(map(x -> real(x[1]), sol.u), map(x -> imag(x[1]), sol.u), zcolor=sol.t, xlims=(-8, 8), ylims=(-8, 8), colormap=:viridis, markerstrokewidth=0, markersize=0.8, label="N_A=$(N_A)", colorbar=false)
    push!(plot_list, p)
end

combined_plot = plot(plot_list..., layout=(2, 3), size=(1800, 1200))
savefig(combined_plot, "traj_buildup_NA_combined_plot.png")

plot_list = []

for ω_tilde in ω_tildes
    for P in Ps
        sol = run_dynamics_with_hetero(P, ω_tilde)
        p = scatter(map(x -> real(x[1]), sol.u), map(x -> imag(x[1]), sol.u), zcolor=sol.t, xlims=(-8, 8), ylims=(-8, 8), colormap=:viridis, markerstrokewidth=0, markersize=0.8, label="P=$(P), ω=$(ω_tilde)", colorbar=false)
        push!(plot_list, p)
    end
end

combined_plot = plot(plot_list..., layout=(length(Ps), length(ω_tildes)), size=(2600, 1600))
xaxis!(combined_plot, L"Re[\lambda]")
yaxis!(combined_plot, L"Im[\lambda]")
savefig(combined_plot, "traj_buildup_combined_plot.png")