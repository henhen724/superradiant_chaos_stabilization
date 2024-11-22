using DifferentialEquations, Plots

# Define the parameters
η_x = 0.1
η_y = 0.1

# Define the drift and diffusion functions
function drift!(du, u, p, t)
    du[1] = u[1]^2 - u[2]^2
    du[2] = -2 * u[2]
end

# Define the vector field for the drift function
function vector_field(x, y)
    u = [x, y]
    du = similar(u)
    drift!(du, u, nothing, 0.0)
    return du
end

# Generate a grid of points
x = y = range(-1.0, 1.0, length=20)
x, y = [range(-1.0, 1.0, length=20) for _ in 1:2]
u = [0.1 * vector_field(xi, yi) for xi in x, yi in y]

# Extract the components of the vector field
u1 = [u[i, j][1] for i in 1:length(x), j in 1:length(y)]
u2 = [u[i, j][2] for i in 1:length(x), j in 1:length(y)]

# Plot the vector field
quiver(x, y, quiver=(u1, u2), title="Drift Vector Field", xlabel="x", ylabel="y")

function diffusion!(du, u, p, t)
    du[1, 1] = η_x
    du[2, 2] = η_y
end

# Initial conditions
u0 = [0.01, 0.01]

# Time span
tspan = (0.0, 10.0)

# Define the SDE problem
prob = SDEProblem(drift!, diffusion!, u0, tspan)

# Solve the SDE problem
sol = solve(prob)

# Plot the solution
using Plots
plot(sol, vars=(1, 2), title="SDE Solution", xlabel="Time", ylabel="Values")