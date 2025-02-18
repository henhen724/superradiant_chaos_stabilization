using QuantumOptics, ForwardDiff, Plots, Optim

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

# Function to create a displacement operator
function displacement(alpha, a_op)
    return exp(-im * (conj(alpha) * a_op + alpha * dagger(a_op)))
end

# Function to calculate the finite difference gradient of the displacement operator
function finite_difference_gradient(alpha, a_op, h=1e-5)
    D_op = displacement(alpha, a_op)
    D_op_plus = displacement(alpha + h, a_op)
    D_op_minus = displacement(alpha - h, a_op)
    grad_real = (D_op_plus - D_op_minus) / (2 * h)

    D_op_plus_im = displacement(alpha + im * h, a_op)
    D_op_minus_im = displacement(alpha - im * h, a_op)
    grad_imag = (D_op_plus_im - D_op_minus_im) / (2 * h)

    return grad_real, grad_imag
end

# Example usage of the finite difference gradient function
alpha = 0.1 + 0.1im
grad_real, grad_imag = finite_difference_gradient(alpha, a_op)
println("Finite difference gradient (real part):")
println(grad_real)
println("Finite difference gradient (imaginary part):")
println(grad_imag)

# Calculate the displacement operator
D_op = displacement(alpha, a_op)

# Print the displacement operator
println("Displacement operator D(α):")
println(D_op)

# Define a vector of 2-tuples of complex numbers
alphas = [0.1]
betas = [1.1]

function op_evolution(op_0, alphas, betas)
    op_t = [op_0]
    for i in eachindex(alphas)
        U_q = displacement(alphas[i], σ_minus / 2)
        U_qc = displacement(betas[i], dagger(a_op) * σ_minus / 2)
        new_op = U_q * U_qc * op_t[end] * dagger(U_qc) * dagger(U_q)
        push!(op_t, new_op)
    end
    return op_t
end

function op_evolution_rev(op_0, alphas, betas)
    op_t = [op_0]
    for i in eachindex(alphas)
        U_q = displacement(alphas[i], σ_minus / 2)
        U_qc = displacement(betas[i], dagger(a_op) * σ_minus / 2)
        new_op = dagger(U_qc) * dagger(U_q) * op_t[end] * U_q * U_qc
        push!(op_t, new_op)
    end
    return reverse(op_t)
end

function calculate_grape_gradient(ρ_t, σ_t, alphas, betas)
    alpha_grad = zeros(eltype(alphas), length(alphas))
    beta_grad = zeros(eltype(betas), length(betas))
    for i in eachindex(alphas)
        U_qc = displacement(betas[i], dagger(a_op) * σ_minus / 2)
        U_q = displacement(alphas[i], σ_minus / 2)

        # Calculate the derivative of the unitary with respect to the real part of the displacement parameter
        dU_qc_dbeta = (-im * betas[i] * a_op * dagger(σ_minus) / 2) * displacement(betas[i], dagger(a_op) * σ_minus / 2)
        dU_q_dalpha = (-im * dagger(σ_minus) / 2) * displacement(alphas[i], σ_minus / 2)

        commutator = dU_qc_dbeta * ρ_t[i] - ρ_t[i] * dU_qc_dbeta
        beta_grad[i] = tr(dagger(σ_t[i+1]) * commutator)
        mid_ρ = U_qc * ρ_t[i] * dagger(U_qc)
        commutator = dU_q_dalpha * mid_ρ - mid_ρ * dU_q_dalpha
        mid_σ = dagger(U_q) * σ_t[i+1] * U_qc
        alpha_grad[i] = tr(dagger(mid_σ) * commutator)
    end
    return alpha_grad, beta_grad
end

# Example usage of the apply_control function
initial_state = tensor(fockstate(b_fock, 0), spindown(b_spin))
ρ_t = op_evolution(dm(initial_state), alphas, betas)
target_state = tensor(fockstate(b_fock, 1), spindown(b_spin))
σ_t = op_evolution_rev(dm(target_state), alphas, betas)

alpha_grad, beta_grad = calculate_grape_gradient(ρ_t, σ_t, alphas, betas)

# Define the objective function
function objective(x)
    alphas = x[1:length(x)÷2]
    betas = x[length(x)÷2+1:end]

    ρ_t = op_evolution(dm(initial_state), alphas, betas)
    final_state = ρ_t[end]

    fidelity = abs(tr(dagger(dm(target_state)) * final_state))^2
    return 1 - fidelity
end

# Define the gradient function
function gradient!(G, x)
    alphas = x[1:length(x)÷2]
    betas = x[length(x)÷2+1:end]

    ρ_t = op_evolution(dm(initial_state), alphas, betas)
    σ_t = op_evolution_rev(dm(target_state), alphas, betas)

    alpha_grad, beta_grad = calculate_grape_gradient(ρ_t, σ_t, alphas, betas)
    G[1:length(x)÷2] = alpha_grad
    G[length(x)÷2+1:end] = beta_grad
end

# Optimize using LBFGS
x0 = vcat(alphas, betas)
result = optimize(objective, gradient!, x0, LBFGS(), Optim.Options(g_tol=1e-5, show_trace=true, iterations=1000))

# Extract optimized alphas and betas
optimized_alphas = result.minimizer[1:length(x0)÷2]
optimized_betas = result.minimizer[length(x0)÷2+1:end]

# Plot the cost over time
costs = [objective(vcat(alphas, betas)) for _ in 1:1000]
plot(1:1000, costs, xlabel="Iteration", ylabel="Cost", title="Cost over Time")

ρ_t_boson = [ptrace(ρ, 2) for ρ in ρ_t]
ρ_t_spin = [ptrace(ρ, 1) for ρ in ρ_t]

function plot_wigner(ρ_t::Vector{<:Operator}; args...)
    wigner_max_and_min = [
        let W = wigner(ρ, range(-5, stop=5, length=100), range(-5, stop=5, length=100))
            (maximum(W), minimum(W))
        end for ρ in ρ_t
    ]
    wigner_max = maximum([max_t for (max_t, min_t) in wigner_max_and_min])
    wigner_min = minimum([min_t for (max_t, min_t) in wigner_max_and_min])

    anim = @animate for ρ in ρ_t
        plot_wigner(ρ; clims=(wigner_min, wigner_max), args...)
    end

    gif(anim, "wigner_animation.gif", fps=10)
end

function plot_wigner(ρ::Operator; c=:viridis, xlabel="x", ylabel="p", title="Wigner Function", kwargs...)
    x = y = range(-5, stop=5, length=100)
    W = wigner(ρ, x, y)

    heatmap(x, y, W; c=c, xlabel=xlabel, ylabel=ylabel, title=title, kwargs...)
end

plot_wigner(ρ_t_boson)