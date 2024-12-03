using LinearAlgebra
using NLsolve
using FiniteDiff

# Define your system of equations
function f(x)
    return [
        x[1] * (1 - x[2]),
        x[2] * (x[1] - 1)
    ]
end

function main()
    # Initial guess
    initial_guess = [0.5, 0.5]

    # Numerical root finding
    result = nlsolve(f, initial_guess)
    fixed_point = result.zero

    println("Fixed Point: ", fixed_point)

    # Compute the numerical Jacobian at the fixed point
    function jacobian(f, x)
        J = FiniteDiff.finite_difference_jacobian(f, x)
        return J
    end

    jacobian_matrix = jacobian(f, fixed_point)
    println("Jacobian Matrix at fixed point: ", jacobian_matrix)

    # Compute eigenvalues
    eigenvalues = eigen(jacobian_matrix).values
    println("Eigenvalues: ", eigenvalues)

    # Determine stability
    stable = all(real(eigenvalues) .< 0)
    unstable = any(real(eigenvalues) .> 0)

    println("Is the fixed point stable? ", stable)
    println("Is the fixed point unstable? ", unstable)
end

main()