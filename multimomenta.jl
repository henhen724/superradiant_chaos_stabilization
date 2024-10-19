function to_1d_index(i::Int, j::Int, ncols::Int)::Int
    return (i - 1) * ncols + j
end

function to_2d_indices(index::Int, ncols::Int)::Tuple{Int,Int}
    i = div(index - 1, ncols) + 1
    j = mod(index - 1, ncols) + 1
    return (i, j)
end

function safe_index(vec::Vector{T}, idx::Int) where {T}
    return (1 <= idx <= length(vec)) ? vec[idx] : zero(T)
end

function safe_index_2D(vec::Vector{T}, i::Int, j::Int, nrows::Int, ncols::Int) where {T}
    if 1 <= i <= nrows && 1 <= j <= ncols
        return vec[1+to_1d_index(i, j, ncols)]
    else
        return zero(T)
    end
end

using DifferentialEquations

momentum_cutoff_long = 21
momentum_cutoff_trans = 21
@assert momentum_cutoff_long % 2 == 1 & momentum_cutoff_trans % 2 == 1
vec_dim = 1 + momentum_cutoff_long * momentum_cutoff_trans
u0 = zeros(ComplexF64, vec_dim)

noise_prototype = zeros(ComplexF64, (vec_dim, 2)) #hetrodyne
noise = RealWienerProcess!(0.0, zeros(2), save_everystep=false)
# noise_prototype = zeros(Float64, (vec_dim, 1)) #homodyne

N_A = 10^5 #unitless
u0[2] = sqrt(N_A)
κ = 8.1 #MHz
E_0 = 40.0 #MHz
ω_r = 0.05 #MHz


ω_c = 80.0 #MHz
P = 1.0 #MHz



function sde_drift(du, u, p, t)
    long_sum_all = 0
    checker_board_all = 0
    for i in 1:momentum_cutoff_long
        for j in 1:momentum_cutoff_trans
            mom_indx = to_1d_index(i, j, momentum_cutoff_trans)
            n = i - (momentum_cutoff_long + 1) // 2
            m = j - (momentum_cutoff_trans + 1) // 2
            trans_sum = safe_index_2D(u, i, j + 2, momentum_cutoff_long, momentum_cutoff_trans) + safe_index_2D(u, i, j - 2, momentum_cutoff_long, momentum_cutoff_trans)
            long_sum = safe_index_2D(u, i + 2, j, momentum_cutoff_long, momentum_cutoff_trans) + safe_index_2D(u, i - 2, j, momentum_cutoff_long, momentum_cutoff_trans)
            long_sum_all += conj(long_sum) * u[1+mom_indx]
            checker_board = safe_index_2D(u, i + 1, j + 1, momentum_cutoff_long, momentum_cutoff_trans) +
                            safe_index_2D(u, i + 1, j - 1, momentum_cutoff_long, momentum_cutoff_trans) +
                            safe_index_2D(u, i - 1, j + 1, momentum_cutoff_long, momentum_cutoff_trans) +
                            safe_index_2D(u, i - 1, j - 1, momentum_cutoff_long, momentum_cutoff_trans)
            checker_board_all += conj(checker_board) * u[1+mom_indx]
            du[1+mom_indx] = -im * ω_r * ((n^2 + m^2) * u[1+mom_indx] - conj(u[1]) * u[1] * long_sum - P * (u[1] + conj(u[1])) * checker_board - P^2 * trans_sum)
        end
    end
    du[1] = -(κ + im * ω_c) * u[1] + im * E_0 * u[1] * long_sum_all + im * E_0 * P * checker_board_all
end

function sde_diffusion(du, u, p, t)
    du[1, 1] = sqrt(κ / 4)
    du[1, 2] = im * sqrt(κ / 4) # comment out for homodyne
end


tspan = (0.0, 1.0)

prob = SDEProblem(sde_drift, sde_diffusion, u0, tspan, noise_rate_prototype=noise_prototype, noise=noise)
sol = solve(prob, RKMilGeneral(; ii_approx=IICommutative()))

using Plots
plot(sol.t, map(x -> real(x[1]), sol.u))