using DifferentialEquations, Plots, NonlinearSolve

# Define the 2D ODE system
function odesystem!(du, u, p, t)
    du[1] = -im * (abs(u[1])^2 + abs(u[2])^2 - 4) * u[1]
    du[2] = -(abs(u[1])^2 + abs(u[2])^2 - 4) * u[2]
end

# Initial conditions
u0 = ComplexF64[1.0, 0.1]

# Time span
tspan = (0.0, 10.0)

# Define the problem
prob = ODEProblem(odesystem!, u0, tspan)

# Solve the problem
sol = solve(prob)

# Plot the solution
plot(map(x -> real(x[1]), sol.u), map(x -> real(x[2]), sol.u), xlabel="Time", ylabel="u", title="Solution of the ODE System")

NLprob = SteadyStateProblem(odesystem!, u0)
NLsol = solve(NLprob; abstol=1e-5, reltol=1e-5)