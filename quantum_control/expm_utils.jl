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


function expm_der_svd(A, B; N_cutoff=nothing)
    U, s, V = svd(A)

    B_tilde = V * B * U

    result = Matrix([s[i] == s[j] ? B_tilde[i, j] : (exp(s[j] - s[i]) - 1) / (s[j] - s[i]) * B_tilde[i, j] for i in 1:size(B_tilde, 1), j in 1:size(B_tilde, 2)])

    return exp(A) * result
end

function expm_der_alt(A, B; N_cutoff=4000)
    dt = 1 / N_cutoff
    result = zero(A)
    for i in 0:N_cutoff
        t = i * dt
        result += exp(-t * A) * B * exp(t * A) * dt
    end

    return exp(A) * result
end

# Function to create a displacement operator
function displacement(params, a_op)
    r = sqrt(params[1]^2 + params[2]^2)
    new_r = r % 4π
    alpha = new_r / r * complex(params[1], params[2])
    # alpha = complex(params[1], params[2])
    return exp(-im * (conj(alpha) * a_op + alpha * dagger(a_op)))
end

function displacement(alpha::Complex, a_op)
    return exp(-im * (conj(alpha) * a_op + alpha * dagger(a_op)))
end

function displacement_dparam(params, a_op)
    r = sqrt(params[1]^2 + params[2]^2)
    new_r = r % 4π
    alpha = new_r / r * complex(params[1], params[2])
    alpha = complex(params[1], params[2])
    return [-im * expm_der(-im * (conj(alpha) * a_op + alpha * dagger(a_op)), a_op + dagger(a_op)), expm_der(-im * (conj(alpha) * a_op + alpha * dagger(a_op)), -a_op + dagger(a_op))]
end