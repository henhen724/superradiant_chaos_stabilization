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

using DifferentialEquations

longmax = 1
transmax = 1

N_A = 10^5 #unitless
κ = 8.1 #MHz
E_0 = 40.0 #MHz
ω_r = 0.05 #MHz


ω_c = 1.170 * E_0 #MHz
P = 0.617 #MHz 


vec_dim = 1 + (2 * longmax + 1) * (2 * transmax + 1)
u0 = (1 - im) / (2 * sqrt(2 * N_A)) * ones(ComplexF64, vec_dim)
u0[2+to_1d_index(0, 0, transmax, longmax)] = 1.0
u0[2+to_1d_index(1, 0, transmax, longmax)] = 0.0
u0[2+to_1d_index(0, 1, transmax, longmax)] = 0.0
u0[2+to_1d_index(-1, 0, transmax, longmax)] = 0.0
u0[2+to_1d_index(0, -1, transmax, longmax)] = 0.0

noise_prototype = zeros(ComplexF64, (vec_dim, 2)) #hetrodyne
noise = RealWienerProcess!(0.0, zeros(2), save_everystep=false)
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
    du[1, 1] = sqrt(κ / 4)
    du[1, 2] = im * sqrt(κ / 4) # comment out for homodyne
end


tspan = (0.0, 500000.0)

probSDE = SDEProblem(sde_drift, sde_diffusion, u0, tspan, noise_rate_prototype=noise_prototype, noise=noise)
probODE = ODEProblem(sde_drift, u0, tspan)
# sol = solve(probSDE, RKMilGeneral(; ii_approx=IICommutative()); adaptive=false, dt=2^(-15))
sol = solve(probODE, Tsit5(); reltol=10^-6, abstol=10^-5, dt=10^(-3), maxiters=10^10)
using Plots
plot(sol.t, map(x -> real(x[1]), sol.u))
plot(map(x -> real(x[1]), sol.u)[end-200:end], map(x -> imag(x[1]), sol.u)[end-200:end])