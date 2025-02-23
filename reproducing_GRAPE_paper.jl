using QuantumOptics, Plots, Optim, FiniteDiff, Random, LaTeXStrings
include("plotting_lib.jl")
include("quantum_control/circuit_evolution.jl")

# Define the cutoff for the Fock space
Nfock = 4

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

function expm(A; N_cutoff=20)
    return sum([A^n / factorial(n) for n in 0:N_cutoff])
end

function expm_der(A, B; N_cutoff=20)
    curr_op = B
    factorial = 1
    result = B

    for n in 1:N_cutoff
        factorial *= n + 1
        curr_op = A * curr_op - curr_op * A
        result += curr_op / factorial
    end

    return result * exp(A)
end

# Function to create a displacement operator
function displacement(params, a_op)
    alpha = complex(params[1], params[2])
    return exp(-im * (conj(alpha) * a_op + alpha * dagger(a_op)))
end

function displacement(alpha::Complex, a_op)
    return exp(-im * (conj(alpha) * a_op + alpha * dagger(a_op)))
end

function displacement_dparam(params, a_op)
    alpha = complex(params[1], params[2])
    return [-im * expm_der(-im * (conj(alpha) * a_op + alpha * dagger(a_op)), a_op + dagger(a_op)), expm_der(-im * (conj(alpha) * a_op + alpha * dagger(a_op)), -a_op + dagger(a_op))]
end


init_state = tensor(fockstate(b_fock, 0), spindown(b_spin))
targ_state = tensor(fockstate(b_fock, 3), spindown(b_spin))

# Check that op_evolution and the evolve QuantumCircuit result in the same thing
U_q = Gate(params -> displacement(params, σ_minus / 2), params -> displacement_dparam(params, σ_minus / 2), 2, L"$U_q$")
U_qc = Gate(params -> displacement(params, dagger(a_op) * σ_minus / 2), params -> displacement_dparam(params, dagger(a_op) * σ_minus / 2), 2, L"$U_{qc}$")

UnitCircuit = QuantumCircuit(repeat([U_q, U_qc], outer=3))

draw_circuit(UnitCircuit)

fidelity = Cost((params, ρ) -> 0, ρ -> 1 - real(tr(dagger(dm(targ_state)) * ρ)), (params, ρ) -> 0 * ρ, ρ -> -dm(targ_state), (params, ρ) -> real.(0 * params))

params, result = GRAPE(dm(init_state), UnitCircuit, fidelity)

draw_circuit(UnitCircuit; params=params)

# Calculate the evolution of the initial state using the optimized parameters
ρ_t, _ = evolve(dm(init_state), UnitCircuit, params)

ρ_t_boson = [ptrace(ρ, 2) for ρ in ρ_t]
ρ_t_spin = [ptrace(ρ, 1) for ρ in ρ_t]

plot_wigner(ρ_t_boson[end])

# Plot the trace of the cost function over time
cost_trace = [fidelity(params, ρ) for ρ in ρ_t]

plot(cost_trace, title="Cost Function Trace Over Time", xlabel="Time Step", ylabel="Cost")