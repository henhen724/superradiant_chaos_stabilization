using QuantumOptics, Plots, Optim, FiniteDiff, Random, Test

# Function to tensor an operator with identity operators in all other bases
function mb(op, index, bases)
    ops = [index == i ? op : identityoperator(bases[i]) for i in eachindex(bases)]
    return tensor(ops...)
end

b_spin = SpinBasis(1 // 2)

# Test displacement_dparam function
function test_displacement_dparam(params)
    σ_minus = sigmam(b_spin)
    dparams = displacement_dparam(params, σ_minus)
    alpha = complex(params[1], params[2])
    disp_op = displacement(params, σ_minus)

    # Check the derivatives with respect to the real part of alpha
    delta = 1e-5
    params_real = [params[1] + delta, params[2]]
    disp_op_real = displacement(params_real, σ_minus)
    numerical_derivative_real = (disp_op_real - disp_op) / delta
    @test isapprox(dense(dparams[1]), dense(numerical_derivative_real), atol=1e-2)

    # Check the derivatives with respect to the imaginary part of alpha
    params_imag = [params[1], params[2] + delta]
    disp_op_imag = displacement(params_imag, σ_minus)
    numerical_derivative_imag = (disp_op_imag - disp_op) / delta
    @test isapprox(dense(dparams[2]), dense(numerical_derivative_imag), atol=1e-2)
end

# Run the test
test_displacement_dparam([1.0, 0.3])
test_displacement_dparam([π, 0.3])
test_displacement_dparam([π, 0.0])
test_displacement_dparam([0.0, π])