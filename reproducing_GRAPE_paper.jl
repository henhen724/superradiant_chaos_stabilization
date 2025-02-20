using QuantumOptics, Plots, Optim, FiniteDiff, Random
include("plotting_lib.jl")
include("quantum_control/evolution.jl")

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
    ops = [index == i ? op : identityoperator(bases[i]) for i in 1:length(bases)]
    return tensor(ops...)
end

# Example usage of the mb function
a_op = mb(destroy(b_fock), 1, [b_fock, b_spin])
σ_minus = mb(sigmam(b_spin), 2, [b_fock, b_spin])

function adj(A)
    function commutator(X)
        return A * X - X * A
    end
    return LinearOperator(commutator, size(A))
end

function expm_der(A, B; N_cutoff=20)
    expA = exp(A)
    curr_op = B
    factorial = 1
    result = B

    for n in 1:N_cutoff
        factorial *= n + 1
        curr_op = A * curr_op - curr_op * A
        result += curr_op / factorial
    end

    return expA * result
end

# Function to create a displacement operator
function displacement(alpha, a_op)
    return exp(-im * (conj(alpha) * a_op + alpha * dagger(a_op)))
end

U_q = Gate((params) -> displacement(params[1] + im * params[2], σ_minus / 2), 2, "U_q")
U_qc = Gate((params) -> displacement(params[1] + im * params[2], dagger(a_op) * σ_minus / 2), 2, "U_qc")

UnitCircuit = QuantumCircuit(repeat([U_q, U_qc], outer=3))

function op_evolution(op_0, alphas, betas)
    op_t = [op_0]
    for i in eachindex(alphas)
        U_q = displacement(alphas[i], σ_minus / 2)
        U_qc = displacement(betas[i], dagger(a_op) * σ_minus / 2)
        new_op = U_qc * U_q * op_t[end] * dagger(U_q) * dagger(U_qc)
        push!(op_t, new_op)
    end
    return op_t
end

function op_evolution_rev(op_0, alphas, betas)
    op_t = [op_0]
    for i in eachindex(alphas)
        U_q = displacement(alphas[i], σ_minus / 2)
        U_qc = displacement(betas[i], dagger(a_op) * σ_minus / 2)
        new_op = dagger(U_q) * dagger(U_qc) * op_t[end] * U_qc * U_q
        push!(op_t, new_op)
    end
    return reverse(op_t)
end

function optimize_grape(num_steps, initial_state, target_state; seed=nothing)
    Random.seed!(seed)
    # Generate random alphas and betas using random unit normals
    alphas = 2 * π * ComplexF64[randn() + im * randn() for _ in 1:num_steps]
    betas = 2 * π * ComplexF64[randn() + im * randn() for _ in 1:num_steps]

    function calculate_grape_gradient(ρ_t, σ_t, alphas, betas)
        alpha_grad = zeros(eltype(alphas), length(alphas))
        beta_grad = zeros(eltype(betas), length(betas))
        for i in eachindex(alphas)
            function cost(alpha, beta)
                U_qc = displacement(beta, dagger(a_op) * σ_minus / 2)
                U_q = displacement(alpha, σ_minus / 2)

                return 1 - real(tr(dagger(σ_t[i+1]) * U_qc * U_q * ρ_t[i] * dagger(U_q) * dagger(U_qc)))
            end
            grad = FiniteDiff.finite_difference_gradient((x) -> cost(x[1] + im * x[2], x[3] + im * x[4]), [real(alphas[i]), imag(alphas[i]), real(betas[i]), imag(betas[i])])
            alpha_grad[i] = grad[1] + im * grad[2]
            beta_grad[i] = grad[3] + im * grad[4]
        end
        return alpha_grad, beta_grad
    end

    function objective(x)
        alphas = x[1:length(x)÷2]
        betas = x[length(x)÷2+1:end]

        ρ_t = op_evolution(dm(initial_state), alphas, betas)
        final_state = ρ_t[end]

        fidelity = real(tr(dm(target_state) * final_state))
        return 1 - fidelity
    end

    function gradient!(G, x)
        alphas = x[1:length(x)÷2]
        betas = x[length(x)÷2+1:end]

        ρ_t = op_evolution(dm(initial_state), alphas, betas)
        σ_t = op_evolution_rev(dm(target_state), alphas, betas)

        alpha_grad, beta_grad = calculate_grape_gradient(ρ_t, σ_t, alphas, betas)
        G[1:length(x)÷2] = alpha_grad
        G[length(x)÷2+1:end] = beta_grad
    end

    x0 = vcat(alphas, betas)
    result = optimize(objective, gradient!, x0, LBFGS(), Optim.Options(g_tol=1e-5, show_trace=true, iterations=1000))

    optimized_alphas = result.minimizer[1:length(x0)÷2]
    optimized_betas = result.minimizer[length(x0)÷2+1:end]

    return optimized_alphas, optimized_betas
end


init_state = tensor(fockstate(b_fock, 0), spindown(b_spin))
targ_state = tensor(fockstate(b_fock, 3), spindown(b_spin))

# Check that op_evolution and the evolve QuantumCircuit result in the same thing

# Define the parameters for the QuantumCircuit
c_params::Vector{Vector{Float64}} = []
for i in eachindex(optimized_alphas)
    push!(c_params, [real(optimized_alphas[i]), imag(optimized_alphas[i])])
    push!(c_params, [real(optimized_betas[i]), imag(optimized_betas[i])])
end

# Evolve the initial state using the QuantumCircuit
ρ_final, ρ_t, measment_record = evolve(dm(init_state), UnitCircuit, c_params)

# Evolve the initial state using the op_evolution function
op_evolution_result = op_evolution(dm(init_state), optimized_alphas, optimized_betas)

# Compare the final states
println("Final state from QuantumCircuit evolution:")
println(ρ_t[end])

println("Final state from op_evolution:")
println(op_evolution_result[end])

# Check if the final states are approximately equal
println("Are the final states approximately equal?")
println(isapprox(ρ_t[end], op_evolution_result[end], atol=1e-5))

optimized_alphas, optimized_betas = optimize_grape(3, init_state, targ_state)

ρ_t = op_evolution(dm(init_state), optimized_alphas, optimized_betas)

ρ_t_boson = [ptrace(ρ, 2) for ρ in ρ_t]
ρ_t_spin = [ptrace(ρ, 1) for ρ in ρ_t]

plot_wigner(ρ_t_boson[end])