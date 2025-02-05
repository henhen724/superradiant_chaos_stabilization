using DifferentialEquations, Plots, BoundaryValueDiffEq

# Define the function f
f(x) = exp(-x^2)
gradf(x) = -2x * exp(-x^2)
gradgradf(x) = -2 * exp(-x^2) + 4x^2 * exp(-x^2)

# Define the ODE
function ode!(du, u, p, t)
    x, λ = u
    du[1] = gradf(x) - λ / 2
    du[2] = (λ - 2 * gradf(x)) * gradgradf(x)
end

u0 = [1.0, -0.1]
tspan = (0.0, 20.0)

function bc1!(residual, sol, p, t)
    residual[1] = sol(0.0)[1] - 1.0 # the solution at the middle of the time span should be -pi/2
    residual[2] = sol(tspan[2])[2]
end

# Initial conditions


# Solve the system of equations
prob = BVProblem(ode!, bc1!, u0, tspan)
sol = solve(prob, MIRK4(), dt=0.01)

# Plot the results
plot(sol, idxs=(0, 1), label="x(t)", xlabel="t", ylabel="x")
plot!(sol.t, -map(x -> x[2], sol.u) / 2, label="u(t) = -λ(t)/2", xlabel="t", ylabel="u")
plot(x -> -(2x * exp(-x^2))^2)