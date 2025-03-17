using QuantumOptics, Plots, Optim, FiniteDiff, Random, LaTeXStrings
include("plotting_lib.jl")
include("quantum_control/circuit_evolution.jl")
include("quantum_control/expm_utils.jl")

# Define the cutoff for the Fock space
Nfock = 30

# Create the Fock space basis
b_fock = FockBasis(Nfock)
b_spin = SpinBasis(1 // 2)

# Create a composite Hilbert space
bases = [b_fock, b_spin]
full_basis = tensor(bases...)

# Function to tensor an operator with identity operators in all other bases
function mb(op, index, bases)
    ops = [index == i ? op : identityoperator(bases[i]) for i in eachindex(bases)]
    return tensor(ops...)
end

# Example usage of the mb function
a_op = mb(destroy(b_fock), 1, [b_fock, b_spin])
σ_minus = mb(sigmam(b_spin), 2, [b_fock, b_spin])


init_state = tensor(fockstate(b_fock, 0), spindown(b_spin))
targ_state = tensor(fockstate(b_fock, 3), spindown(b_spin))

# Check that op_evolution and the evolve QuantumCircuit result in the same thing
U_q = Gate(params -> displacement(params, σ_minus / 2), params -> displacement_dparam(params, σ_minus / 2), nothing, 2, L"$U_q$")
U_qc = Gate(params -> displacement(params, dagger(a_op) * σ_minus / 2), params -> displacement_dparam(params, dagger(a_op) * σ_minus / 2), nothing, 2, L"$U_{qc}$")

UnitCircuit = QuantumCircuit(repeat([U_q, U_qc], outer=3))
fidelity_3 = Cost((params, ρ) -> 0, ρ -> 1 - real(tr(dagger(dm(targ_state)) * ρ)), (params, ρ) -> 0 * ρ, ρ -> -dm(targ_state), (params, ρ) -> real.(0 * params))
params, result = GRAPE(dm(init_state), UnitCircuit, fidelity_3)

# Set the resolution of the figure to 300 dpi
default(dpi=300)
draw_circuit(UnitCircuit; params=params)
savefig("circuit_diagram.png")

# Calculate the evolution of the initial state using the optimized parameters
ρ_t, _ = evolve(dm(init_state), UnitCircuit, params)

ρ_t_boson = [ptrace(ρ, 2) for ρ in ρ_t]
ρ_t_spin = [ptrace(ρ, 1) for ρ in ρ_t]

anim = plot_wigner(ρ_t_boson)
# Plot the photon number over the course of the circuit
photon_numbers = [real(expect(number(b_fock), ρ)) for ρ in ρ_t_boson]

plot(photon_numbers, label="Photon Number", xlabel="Time Step", ylabel="Photon Number", legend=:topright)
savefig("photon_number_plot.png")

# Calculate the excited state population of the atom over the course of the circuit
excited_populations = [real(expect(sigmaz(b_spin), ρ)) for ρ in ρ_t_spin]

# Plot the excited state population
plot(excited_populations, label=L"S_z", xlabel="Time Step", ylabel=L"\langle S_z \rangle", legend=:topright)
savefig("excited_population_plot.png")

targ_state = tensor((coherentstate(b_fock, 3.0) + coherentstate(b_fock, -3.0) + coherentstate(b_fock, 3.0im) + coherentstate(b_fock, -3.0im)), spindown(b_spin))
targ_state = targ_state / norm(targ_state)
plot_wigner(ptrace(dm(targ_state), 2))
fidelity_cat = Cost((params, ρ) -> 0, ρ -> 1 - real(tr(dagger(dm(targ_state)) * ρ)), (params, ρ) -> 0 * ρ, ρ -> -dagger(dm(targ_state)), (params, ρ) -> real.(0 * params))

# Prepare the cat state with circuits of length 1 to 50
UnitCircuit = QuantumCircuit(repeat([U_q, U_qc], outer=20))
params, result = GRAPE(dm(init_state), UnitCircuit, fidelity_cat)

# Plot the infidelity over time for the different circuit realizations
infidelities = [map(x -> x.value, result.trace) for result in results]

p = plot()
for (i, infidelity) in enumerate(infidelities)
    plot!(p, 1:length(infidelity), infidelity, color=RGB(i / length(infidelities), 0, 1 - i / length(infidelities)), label=false)
end
xlabel!("Iteration")
ylabel!("Infidelity")

colormap = [RGB(a, 0, 1 - a) for a in range(0, stop=1, length=100)]

# Add a color bar
using CairoMakie
fig = Figure()
Axis(fig[1, 1])
Colorbar(fig[1, 2], limits=(0, 50), colormap=colormap,
    flipaxis=false)
fig


ρ_t, _ = evolve(dm(init_state), UnitCircuit, params)

ρ_t_boson = [ptrace(ρ, 2) for ρ in ρ_t]
ρ_t_spin = [ptrace(ρ, 1) for ρ in ρ_t]

plot_wigner(ρ_t_boson[end])

iterations = 1:length(result.trace)
g_norms = map(x -> x.g_norm, result.trace)
values = map(x -> x.value, result.trace)

plot(iterations, g_norms, label="g_norm", xlabel="Iteration", ylabel="Value", legend=:topright)
plot(iterations, values, label="Value")