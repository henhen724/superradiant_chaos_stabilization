using QuantumOptics, DifferentialEquations, Plots

# Define parameters
N = 2  # Number of spins

# Define spin operators
basis = SpinBasis(1 // 2)
idOp = identityoperator(basis)
S_z = sigmaz(basis)
S_x = sigmax(basis)
S_y = sigmay(basis)

H_list = [(S_z^2).data, S_x.data, S_z.data]
controls = zeros(Float64, length(H_list))
control_cost = [1.0, 1.0, 1.0]

# Initial wave functions
ψ0 = exp(-im * (0.3 * π / 4) * S_x) * spindown(basis)
ψ0_d = ψ0.data
χ0 = spindown(basis)
u0 = vcat(ψ0.data, χ0.data)#zeros(ComplexF64, length(basis)))
target = (spinup(basis) + spindown(basis)) / sqrt(2)
target_d = target.data

HS_dim = length(basis)
# Define the Schrödinger equation
function quad_pontryagin!(dv, v, u, t)
    dψ, dχ = view(dv, 1:HS_dim), view(dv, HS_dim+1:2*HS_dim)
    ψ, χ = v[1:HS_dim], v[HS_dim+1:2*HS_dim]

    H_tot = 0.0 * H_list[1]
    for (i, H) in enumerate(H_list)
        controls[i] = imag(conj(transpose(χ)) * H * ψ)
        H_tot += controls[i] * H
    end
    dψ .= -im * H_tot * ψ
    dχ .= -im * H_tot * χ
end

T_final = 1.0
tspan = (0.0, T_final)

function bc!(residual, u, p, t)
    residual[1:HS_dim] = u[begin][1:HS_dim] - ψ0_d
    residual[HS_dim+1:end] = u[end][HS_dim+1:end] + target_d * (1 - conj(transpose(target_d)) * u[end][1:HS_dim])
end

# Time span

prob = BVProblem(quad_pontryagin!, bc!, u0, tspan)

function normalize_callback!(integrator)
    u = integrator.u
    ψ = view(u, 1:HS_dim)
    norm_factor = sqrt(sum(abs2, ψ))
    ψ .= ψ / norm_factor
end

cb = CallbackSet(DiscreteCallback((u, t, integrator) -> true, normalize_callback!))

sol = solve(prob, MIRK4(), dt=0.01, callback=cb)

# # Define the ODE problem
# prob = ODEProblem(quad_pontryagin!, u0, tspan)

# # Create the integrator
# integrator = init(prob, Tsit5())

# # Integrate the ODE
# while integrator.t < tspan[2]
#     step!(integrator)
# end

# # Extract the solution
# sol = integrator.sol


print("Initial Fidelity: ", abs(conj(transpose(target_d)) * sol.u[begin][1:HS_dim])^2)
print("Final Fidelity: ", abs(conj(transpose(target_d)) * sol.u[end][1:HS_dim])^2)

# Convert the solution into QuantumOptics states
states = [Ket(basis, view(sol.u[i], 1:HS_dim)) for i in 1:length(sol.u)]

# Plot the normalization of the wavefunction over time
norms = [sqrt(sum(abs2, view(sol.u[i], 1:HS_dim))) for i in 1:length(sol.u)]
plot(sol.t, norms, xlabel="Time", ylabel="Normalization", label="|ψ|", title="Normalization of the Wavefunction")

states = states ./ expect(idOp, states)

# Plot the fidelity against the target state over time
plot(sol.t, real.(expect(dm(target), states)), xlabel="Time", ylabel="Fidelity", label="Fidelity", title="Fidelity of Cat State over Time")



plot(sol.t, real.(expect(S_x, states)), xlabel="Time", ylabel="Expectation value", label="S_x", title="Evolution under control")
plot!(sol.t, real.(expect(S_y, states)), xlabel="Time", ylabel="Expectation value", label="S_y")
plot!(sol.t, real.(expect(S_z, states)), xlabel="Time", ylabel="Expectation value", label="S_z")