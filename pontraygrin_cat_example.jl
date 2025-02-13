using QuantumOptics, DifferentialEquations, Plots, LinearAlgebra, LaTeXStrings, Optim

# Define parameters
N = 1  # Number of spins
dt = 0.01  # Time step
T_final = 1.0  # Final time
num_steps = Int(T_final / dt) + 1  # Number of time steps
learning_rate = 0.001  # Learning rate for gradient ascent

# Define spin operators
basis = SpinBasis(N // 2)
idOp = identityoperator(basis)
S_z = sigmaz(basis)
S_x = sigmax(basis)
S_y = sigmay(basis)

ω_z = 0.012 # MHz
H_0 = ω_z * S_z
H_list = [(S_x^2).data, S_x.data, (S_z^2).data, S_z.data]
control_labels = [L"S_x^2", L"S_x", L"S_z^2", L"S_z"]
controls = cumsum(randn(Float64, length(H_list), num_steps) * sqrt(dt), dims=2)
control_cost = [1.0, 1.0, 1.0, 1.0]

# Initial wave functions
ψ0 = exp(-im * (0.01 * π / 4) * S_x) * spindown(basis)
ψ0_d = ψ0.data

target = (spinup(basis) + spindown(basis)) / sqrt(2)
target_d = target.data

HS_dim = length(basis)

# Define the Schrödinger equation
function schrodinger!(dψ, ψ, u, t)
    H_tot = 0.0 * H_list[1]
    for (i, H) in enumerate(H_list)
        H_tot += u[i] * H
    end
    dψ .= -im * H_tot * ψ
end

# Define the cost function
function final_cost_function(final_state, target_d)
    fidelity = abs(conj(transpose(target_d)) * final_state)^2
    return 1 - fidelity
end
function total_cost_function(sol, target_d)
    final_cost = final_cost_function(sol.u[end], target_d)
    running_cost = 0.5 * sum(sum(control_cost .* controls .^ 2)) * dt
    return final_cost + running_cost
end


# Compute the gradient of the cost function
function compute_gradient(sol_forward, sol_backward, target_d, controls, dt)
    grad = zeros(size(controls))
    for t in 1:num_steps
        ψ = sol_forward.u[t]
        χ = sol_backward.u[t]
        for (i, H) in enumerate(H_list)
            grad[i, t] = (control_cost[i] * controls[i, t] - imag(conj(transpose(χ)) * H * ψ)) * dt
        end
    end
    return grad
end

function linear_interpolation(x, y, xq)
    result = zeros(length(xq))
    for (i, xi) in enumerate(xq)
        if xi <= x[1]
            result[i] = y[1]
        elseif xi >= x[end]
            result[i] = y[end]
        else
            j = findlast(t -> t <= xi, x)
            t = (xi - x[j]) / (x[j+1] - x[j])
            result[i] = (1 - t) * y[j] + t * y[j+1]
        end
    end
    return result
end

function optimize_controls_duration!(controls, sol_forward, target_d, ψ0_d, dt, T_final)
    # Find timestep with lowest cost
    costs = [final_cost_function(sol_forward(t), target_d) for t in 0:dt:T_final]
    min_cost_index = argmin(costs)

    # Create new time points for interpolation
    old_t = 0:dt:((min_cost_index-1)*dt)
    new_t = 0:dt:T_final

    # Rescale each control channel
    for i in 1:size(controls, 1)
        old_controls = controls[i, 1:min_cost_index]
        controls[i, :] .= ((min_cost_index - 1) * dt / T_final) * linear_interpolation(old_t, old_controls, new_t)
    end

    return controls
end

sol_forward = nothing
sol_backward = nothing
maxiters = 200
cost_tol = 0.01
last_cost = Inf
# GRAPE algorithm
# Define the objective function for Optim
function objective(x)
    controls_reshaped = reshape(x, size(controls))

    # Forward evolution
    prob_forward = ODEProblem(schrodinger!, ψ0_d, (0.0, T_final), controls_reshaped)
    sol_forward = solve(prob_forward, Tsit5(); saveat=0:dt:T_final)

    # Calculate cost
    return total_cost_function(sol_forward, target_d)
end

# Define the gradient function for Optim
function gradient!(G, x)
    controls_reshaped = reshape(x, size(controls))

    # Forward evolution
    prob_forward = ODEProblem(schrodinger!, ψ0_d, (0.0, T_final), controls_reshaped)
    sol_forward = solve(prob_forward, Tsit5(); saveat=0:dt:T_final)

    # Initialize costate
    λT = 2 * target_d * conj(transpose(target_d)) * sol_forward.u[end]

    # Backward evolution
    prob_backward = ODEProblem(schrodinger!, λT, (T_final, 0.0), controls_reshaped)
    sol_backward = solve(prob_backward, Tsit5(); saveat=T_final:-dt:0)

    # Compute gradient
    grad = compute_gradient(sol_forward, sol_backward, target_d, controls_reshaped, dt)
    G[:] = vec(grad)
end

# Optimize using BFGS
# Initialize array to store control histories
control_history = Vector{Matrix{Float64}}()
push!(control_history, copy(controls))

function callback(x)
    # Store current controls
    println(x.metadata["time"])
    push!(control_history, copy(reshape(x, size(controls))))
    return false
end

maxiters = 1000
result = optimize(objective, gradient!, vec(controls), LBFGS(),
    Optim.Options(iterations=maxiters, g_tol=1e-3, show_trace=true, show_every=10, callback=callback))

# Create animation
anim = @animate for i in 1:length(control_history)
    plot(sol_forward.t[1:size(controls, 2)],
        [control_history[i][j, :] for j in 1:length(H_list)]...,
        xaxis="Time",
        yaxis="Control amplitude",
        label=permutedims(control_labels),
        title="Control amplitudes (iteration $i)",
        legend=:outertopright)
end

gif(anim, "control_convergence.gif", fps=10)

# Update controls with optimized result
controls[:] = reshape(Optim.minimizer(result), size(controls))

# Final forward evolution for plotting
prob_forward = ODEProblem(schrodinger!, ψ0_d, (0.0, T_final), controls)
sol_forward = solve(prob_forward, Tsit5(); saveat=0:dt:T_final)

# Final backward evolution for gradient plotting
λT = 2 * target_d * conj(transpose(target_d)) * sol_forward.u[end]
prob_backward = ODEProblem(schrodinger!, λT, (T_final, 0.0), controls)
sol_backward = solve(prob_backward, Tsit5(); saveat=T_final:-dt:0)

print("Initial Fidelity: ", abs(conj(transpose(target_d)) * sol_forward.u[begin])^2)
print("Final Fidelity: ", abs(conj(transpose(target_d)) * sol_forward.u[end])^2)

# Convert the solution into QuantumOptics states
states = [Ket(basis, sol_forward.u[i]) for i in 1:length(sol_forward.u)]

# Plot the normalization of the wavefunction over time
norms = [sqrt(sum(abs2, sol_forward.u[i])) for i in 1:length(sol_forward.u)]
plot(sol_forward.t, norms, xlabel="Time", ylabel="Normalization", label="|ψ|", title="Normalization of the Wavefunction")

states = states ./ expect(idOp, states)

# Plot the fidelity against the target state over time
plot(sol_forward.t, real.(expect(dm(target), states)), xlabel="Time", ylabel="Fidelity", label="Fidelity", title="Fidelity of Cat State over Time")

fig = plot()
plot!(fig, sol_forward.t, real.(expect(S_x, states)), xlabel="Time", ylabel="Expectation value", label=L"S_x", title="Evolution under control")
plot!(fig, sol_forward.t, real.(expect(S_y, states)), xlabel="Time", ylabel="Expectation value", label=L"S_y")
plot!(fig, sol_forward.t, real.(expect(S_z, states)), xlabel="Time", ylabel="Expectation value", label=L"S_z")
display(fig)

fig = plot()
for i in 1:length(H_list)
    plot!(fig, sol_forward.t[1:size(controls, 2)], controls[i, :], xaxis="Time", yaxis="Control amplitude", label="Control $(control_labels[i])", title="Control amplitudes over time")
end
display(fig)


fig = plot()
grad = compute_gradient(sol_forward, sol_backward, target_d, controls, dt)
for i in 1:length(H_list)
    plot!(fig, 1:num_steps, grad[i, :], xaxis="Time step", yaxis="Gradient amplitude", label="Gradient $i", title="Gradient amplitudes over time")
end
display(fig)