using QuantumOptics
using DifferentialEquations

# Define parameters
N = 2  # Number of spins
ω_z = 1.0  # Coefficient for S_z^2 term
ω_x = 1.0  # Coefficient for S_x term

# Define spin operators
basis = SpinBasis(1 // 2)
S_z = sigmaz(basis)
S_x = sigmax(basis)

# Define Hamiltonian
H = ω_z * S_z^2 + ω_x * S_x

# Initial wave functions
ψ1 = spindown(basis)
ψ2 = spindown(basis)

# Define the Schrödinger equation
function schrodinger!(dv, v, H, t)
    dv .= -im * H.data * v
end

# Time span
tspan = (0.0, 10.0)

# Solve for ψ1
prob1 = ODEProblem(schrodinger!, ψ1, tspan, H)
sol1 = solve(prob1)

# Solve for ψ2
prob2 = ODEProblem(schrodinger!, ψ2, tspan, H)
sol2 = solve(prob2)

# Output the solutions
println("Solution for ψ1: ", sol1)
println("Solution for ψ2: ", sol2)