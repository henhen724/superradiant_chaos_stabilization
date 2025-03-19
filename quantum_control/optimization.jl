using Optimization, Zygote
rosenbrock(u, p) = (p[1] - u[1])^2 + p[2] * (u[2] - u[1]^2)^2
u0 = zeros(2)
p = [1.0, 100.0]

optf = OptimizationFunction(rosenbrock, AutoZygote())
prob = OptimizationProblem(optf, u0, p)

sol = solve(prob, Optimization.LBFGS())

using Plots


# Plot the loss over iterations
plot(loss_values, xlabel="Iteration", ylabel="Loss", title="Loss over Iterations", lw=2)