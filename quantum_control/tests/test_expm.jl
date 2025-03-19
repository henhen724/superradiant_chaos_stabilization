using QuantumOptics, Plots, Optim, FiniteDiff, Random, Test
include("../expm_utils.jl")

# Test suite for displacement_dparam function
@testset "Displacement Derivative Tests" begin
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
        disp_op = displacement(params, σ_minus)

        # Check the derivatives with respect to the real part of alpha
        delta = 1e-5
        params_real = [params[1] + delta, params[2]]
        disp_op_real = displacement(params_real, σ_minus)
        numerical_derivative_real = (disp_op_real - disp_op) / delta
        @test isapprox(dense(dparams[1]), dense(numerical_derivative_real), atol=1e-5)

        # Check the derivatives with respect to the imaginary part of alpha
        params_imag = [params[1], params[2] + delta]
        disp_op_imag = displacement(params_imag, σ_minus)
        numerical_derivative_imag = (disp_op_imag - disp_op) / delta
        @test isapprox(dense(dparams[2]), dense(numerical_derivative_imag), atol=1e-5)
    end
    # Test with different parameter sets
    @testset "Test with params [π, 0.3]" begin
        test_displacement_dparam([π, 0.3])
    end

    @testset "Test with params [π, 0.0]" begin
        test_displacement_dparam([π, 0.0])
    end

    @testset "Test with params [0.0, π]" begin
        test_displacement_dparam([0.0, π])
    end

    # Test with 10 random parameter sets
    Random.seed!(1234)  # Set seed for reproducibility
    for i in 1:10
        params = [2π * rand(), 2π * rand()]  # Random parameters in the range [0, 2π]
        @testset "Test with random params $params" begin
            test_displacement_dparam(params)
        end
    end
end

function test_exp_dparam(matrix, dmatrix)
    delta = 1e-5
    # Compute the derivative using expm_dparam
    # expm_derivative = expm_der_phiv(matrix, dmatrix)
    exp_derivative_diag = expm_der_diag(matrix, dmatrix)

    # Compute the numerical derivative
    numerical_matrix = matrix + delta * dmatrix
    numerical_derivative = (exp(numerical_matrix) - exp(matrix)) / delta

    # Check if the computed derivative matches the numerical derivative
    # @test isapprox(dense(expm_derivative), dense(numerical_derivative), atol=1e-5)
    @test isapprox(exp_derivative_diag, numerical_derivative, atol=1e-4)
end

matrix = -im * [0.0 1.0 0.0; 1.0 0.0 1.0; 0.0 1.0 0.0]
dmatrix = -im * [1.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 -1.0]
test_exp_dparam(matrix, dmatrix)

@testset "Matrix Exponential Derivative Tests" begin
    # Test suite for expm_dparam function
    function test_exp_dparam(matrix, dmatrix)
        delta = 1e-7
        # Compute the derivative using expm_dparam
        # expm_derivative = expm_der_phiv(matrix, dmatrix)
        exp_derivative_diag = expm_der_diag(matrix, dmatrix)

        # Compute the numerical derivative
        numerical_matrix = matrix + delta * dmatrix
        numerical_derivative = (exp(numerical_matrix) - exp(matrix)) / delta

        # Check if the computed derivative matches the numerical derivative
        # @test isapprox(dense(expm_derivative), dense(numerical_derivative), atol=1e-5)
        @test isapprox(exp_derivative_diag, numerical_derivative, atol=1e-3)
    end

    # Test with a simple 2x2 matrix
    @testset "Test with pauli matrix" begin
        matrix = -im * [0.0 1.0; 1.0 0.0]
        dmatrix = -im * [1.0 0.0; 0.0 -1.0]
        test_exp_dparam(matrix, dmatrix)
    end


    @testset "Test with 2x2 matrix" begin
        matrix = -im * [0.0 1.0; 1.0 0.0]
        dmatrix = -im * [1.0 0.0; 0.0 -1.0]
        test_exp_dparam(matrix, dmatrix)
    end


    # Test with a random 3x3 matrix
    @testset "Test with random 3x3 matrix" begin
        Random.seed!(1234)  # Set seed for reproducibility
        matrix = -im * rand(ComplexF64, 3, 3)
        dmatrix = -im * rand(ComplexF64, 3, 3)
        test_exp_dparam(matrix, dmatrix)
    end

    # Test with a diagonal matrix
    @testset "Test with diagonal matrix" begin
        matrix = -im * Diagonal([1.0, 2.0, 3.0])
        dmatrix = -im * Diagonal([0.1, 0.2, 0.3])
        test_exp_dparam(matrix, dmatrix)
    end

    # Test with 5 random matrices
    for i in 1:5
        @testset "Test with random matrix $i" begin
            Random.seed!(i)  # Set seed for reproducibility
            size = rand(2:5)  # Random size between 2 and 5
            matrix = -im * rand(size, size)
            dmatrix = -im * rand(size, size)
            test_exp_dparam(matrix, dmatrix)
        end
    end
end

# Additional tests for the `ad` function to check consistency with the commutator
@testset "Adjoint (ad) Map Tests" begin
    # Function to test the consistency of `ad` with the commutator
    function test_ad_consistency(A, B)
        # Compute the result of `ad(A, B)`
        ad_result = transpose(reshape(ad(A) * reshape(transpose(B), length(B)), size(B)))

        # Compute the commutator [A, B] = A * B - B * A
        commutator = A * B - B * A

        # Check if `ad(A, B)` matches the commutator
        @test isapprox(ad_result, commutator, atol=1e-5)
    end

    # Test with simple 2x2 matrices
    @testset "Test with 2x2 matrices" begin
        A = [1.0 2.0; 3.0 4.0]
        B = [0.5 1.0; 1.5 2.0]
        test_ad_consistency(A, B)
    end

    # Test with random 3x3 matrices
    @testset "Test with random 3x3 matrices" begin
        Random.seed!(1234)  # Set seed for reproducibility
        A = rand(ComplexF64, 3, 3)
        B = rand(ComplexF64, 3, 3)
        test_ad_consistency(A, B)
    end

    # Test with random 3x3 matrices
    @testset "Test with random 3x3 matrices" begin
        Random.seed!(1234)  # Set seed for reproducibility
        A = rand(3, 3)
        B = rand(3, 3)
        test_ad_consistency(A, B)
    end

    # Test with diagonal matrices
    @testset "Test with diagonal matrices" begin
        A = Diagonal([1.0, 2.0, 3.0])
        B = Diagonal([0.5, 1.5, 2.5])
        test_ad_consistency(A, B)
    end

    # Test with sparse matrices
    @testset "Test with sparse matrices" begin
        A = sparse([1.0 0.0 0.0; 0.0 2.0 0.0; 0.0 0.0 3.0])
        B = sparse([0.5 0.0 0.0; 0.0 1.5 0.0; 0.0 0.0 2.5])
        test_ad_consistency(A, B)
    end

    # Test with Hermitian matrices
    @testset "Test with Hermitian matrices" begin
        A = Hermitian([1.0 2.0; 2.0 3.0])
        B = Hermitian([0.5 1.0; 1.0 1.5])
        test_ad_consistency(A, B)
    end

    # Test with 5 random matrices of varying sizes
    for i in 1:5
        @testset "Test with random matrix pair $i" begin
            Random.seed!(i)  # Set seed for reproducibility
            size = rand(2:5)  # Random size between 2 and 5
            A = rand(size, size)
            B = rand(size, size)
            test_ad_consistency(A, B)
        end
    end
end

b_spin = SpinBasis(1 // 2)
σ_minus = sigmam(b_spin)

matrix = [1.0 2.0; 3.0 4.0]
op = Operator(b_spin, matrix)

params = [0.1, 0.0]
dparams = displacement_dparam(params, σ_minus)

# Define a function to compute disp_op and its derivative
function compute_disp_op_and_derivative(param1, param2)
    params = [param1, param2]
    dparams = displacement_dparam(params, σ_minus)
    disp_op = displacement(params, σ_minus)
    return disp_op.data[1, 1], dparams[1].data[1, 1]
end


# Generate data for plotting
param1_values = 0.0:0.01:10π
param2 = 0.6
disp_op_00 = []
disp_op_00_derivative = []
for param1 in param1_values
    value, derivative = compute_disp_op_and_derivative(param1, param2)
    push!(disp_op_00, value)
    push!(disp_op_00_derivative, derivative)
end
# Plot the results
plot(param1_values, real.(disp_op_00), label="Re(disp_op[1,1])", xlabel="param[1]", ylabel="Value")
plot!(param1_values, real.(disp_op_00_derivative), label="Re(derivative)")

plot(param1_values, imag.(disp_op_00), label="Im(disp_op[1,1])", xlabel="param[1]", ylabel="Value")
plot!(param1_values, imag.(disp_op_00_derivative), label="Im(derivative)")

X = ComplexF64[0.0 1.0; 1.0 0.0]
Z = ComplexF64[0.0 -im; im 0.0]
param1_values = 0.0:0.01:10π
disp_op_00 = []
disp_op_00_derivative = []
der_num = []
for param1 in param1_values
    value = exp(-im * X - im * param1 * Z)
    derivative = expm_der_diag(-im * X - im * param1 * Z, -im * Z)
    delta = 1e-5
    value_delta = exp(-im * X - im * (param1 + delta) * Z)
    push!(disp_op_00, value[1, 2])
    push!(disp_op_00_derivative, derivative[1, 2])
    push!(der_num, ((value_delta[1, 2] - value[1, 2]) / delta))
end

plot(param1_values, real.(disp_op_00), label="Re(disp_op[1,1])", xlabel="param[1]", ylabel="Value")
plot!(param1_values, real.(disp_op_00_derivative), label="Re(derivative)")
plot!(param1_values, real.(der_num), label="Re(numerical derivative)")

plot(param1_values, imag.(disp_op_00), label="Im(disp_op[1,1])", xlabel="param[1]", ylabel="Value")
plot!(param1_values, imag.(disp_op_00_derivative), label="Im(derivative)")
plot!(param1_values, imag.(der_num), label="Im(numerical derivative)")