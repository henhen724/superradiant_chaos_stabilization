using DifferentialEquations, Plots, BoundaryValueDiffEq

# Define the function f
f(x) = exp(-x^2)

# Define the ODE
function ode!(du, u, p, t)
    x, λ = u
    du[1] = -2 * x * exp(-x^2) - λ / 2
    du[2] = (λ - 4 * x * exp(-x^2)) * (-2 * exp(-x^2) + 4 * x^2 * exp(-x^2))
end

function bc1!(residual, u, p, t)
    residual[1] = u[begin][1] - 1.0 # the solution at the middle of the time span should be -pi/2
    residual[2] = u[end][2]
end

# Initial conditions
u0 = [1.0, 0.1]
tspan = (0.0, 20.0)

# Solve the system of equations
prob = BVProblem(ode!, bc1!, u0, tspan)
sol = solve(prob, MIRK4(), dt=0.01)

# Plot the results
plot(sol, vars=(0, 1), label="x(t)", xlabel="t", ylabel="x")
plot!(sol.t, -map(x -> x[2], sol.u) / 2, label="u(t) = -λ(t)/2", xlabel="t", ylabel="u")