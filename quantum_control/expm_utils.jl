using LinearAlgebra

function expm(A; N_cutoff=20)
    return sum([A^n / factorial(n) for n in 0:N_cutoff])
end

function expm_der(A, B; N_cutoff=nothing)
    curr_op = B
    factorial = 1
    result = B

    if N_cutoff === nothing
        norm_A = norm(A)
        N_cutoff = ceil(Int, 5 * exp(1) * norm_A)
        N_cutoff = max(N_cutoff, 15)
    end

    for n in 1:N_cutoff
        if n > 10
            factorial = sqrt(2 * π * (n + 1)) * ((n + 1) / ℯ)^(n + 1)
            # println(norm(curr_op) / factorial)
            if isinf(factorial)
                return result * exp(A)
            end
        else
            factorial *= n + 1
        end
        curr_op = A * curr_op - curr_op * A
        result += curr_op / factorial
    end

    return result * exp(A)
end


function expm_der_diag(A::Operator, B::Operator; tol=1e-6)
    der_m = expm_der_diag(Matrix(A.data), Matrix(B.data); tol=tol)
    return Operator(A.basis_l, A.basis_r, der_m)
end

function expm_der_diag(A::AbstractMatrix, B::AbstractMatrix; tol=1e-6)
    vals, U = eigen(A)
    D = Diagonal([vals...])
    U_dag = conj(transpose(U))

    B_tilde = U_dag * B * U

    result = [abs(vals[i] - vals[j]) < tol ? B_tilde[i, j] : (exp(vals[j] - vals[i]) - 1) / (vals[j] - vals[i]) * B_tilde[i, j] for i in axes(B_tilde, 1), j in axes(B_tilde, 2)]

    return exp(A) * U * result * U_dag
end

function ad(A::AbstractMatrix)
    n = size(A, 1)
    iMat = I(n)
    return kron(A, iMat) - kron(iMat, transpose(A))
end

# function expm_der_phiv(A::Operator, B::Operator)
#     der_m = expm_der_phiv(Matrix(A.data), Matrix(B.data))
#     return Operator(A.basis_l, A.basis_r, der_m)
# end

# function expm_der_phiv(A::AbstractMatrix, B::AbstractMatrix)
#     ad_A = ad(A)
#     term1 = phiv(1.0, ad_A, reshape(transpose(B), length(B)), 1)[:, 2]
#     return exp(A) * transpose(reshape(term1, size(B)))
# end

function expm_der_alt(A, B; N_cutoff=4000)
    dt = 1 / N_cutoff
    result = zero(A)
    for i in 0:N_cutoff
        t = i * dt
        result += exp(t * A) * B * exp(-t * A) * dt
    end

    return exp(A) * result
end

# Function to create a displacement operator
function displacement(params, a_op)
    r = params[1] % 4π
    alpha = r * exp(im * params[2])
    return exp(-im * (conj(alpha) * a_op + alpha * dagger(a_op)))
end

function displacement(alpha::Complex, a_op)
    return exp(-im * (conj(alpha) * a_op + alpha * dagger(a_op)))
end

function displacement_dparam(params, a_op)
    r = params[1] % 4π
    alpha = r * exp(im * params[2])
    return [-im * (exp(-im * params[2]) * a_op + exp(im * params[2]) * dagger(a_op)) * exp(-im * (conj(alpha) * a_op + alpha * dagger(a_op))), expm_der_diag(-im * (conj(alpha) * a_op + alpha * dagger(a_op)), alpha * dagger(a_op) - conj(alpha) * a_op)]
end