using QuantumOptics, LinearAlgebra, SparseArrays

sb = SpinBasis(2 // 2)
Sx = sigmax(sb)
Sy = sigmay(sb)
Sz = sigmaz(sb)
idOp = identityoperator(sb)

function are_ops_linear_independent(ops_list)
    matrix = hcat([reshape(op.data, :) for op in ops_list]...)
    return rank(matrix) == length(ops_list)
end

function normalize_op(op)
    trace_squared = tr(dagger(op) * op)
    return op / sqrt(trace_squared)
end

are_ops_linear_independent([Sx, Sy, Sz, idOp])

function find_dynamical_lie_algebra(ops_list; tol=1e-10)
    if length(ops_list) == 0
        return ops_list
    end

    HS_dim = size(ops_list[1].data, 1)
    generated_ops = ops_list

    # Make all operators traceless
    for i in eachindex(generated_ops)
        generated_ops[i] = generated_ops[i] - (tr(generated_ops[i]) / HS_dim) * idOp
    end

    # Normalize all operators so their trace squared is 1
    for i in eachindex(generated_ops)
        generated_ops[i] = normalize_op(generated_ops[i])
    end

    if !are_ops_linear_independent(generated_ops)
        # Turn the generated_ops into vectors and put those vectors into a matrix
        matrix = hcat([reshape(op.data, :) for op in generated_ops]...)

        # Perform singular value decomposition
        if issparse(matrix)
            matrix = Matrix(matrix)
        end
        U, S, V = svd(matrix)

        # Take all the singular vectors with non-zero singular value
        non_zero_singular_vectors = V[S.>tol, :]

        println(size(non_zero_singular_vectors))

        # Use the projection of vectors onto the original ops list to create a new set of linearly independent vectors
        new_generated_ops = []
        for i in 1:size(non_zero_singular_vectors, 1)
            new_op = sum([non_zero_singular_vectors[i, j] * generated_ops[j] for j in eachindex(generated_ops)])
            push!(new_generated_ops, new_op)
        end

        generated_ops = new_generated_ops
    end

    checked_pair_indices = Set{Tuple{Int,Int}}()

    while length(generated_ops) < HS_dim^2 - 1
        if length(checked_pair_indices) == binomial(length(generated_ops), 2)
            break
        end
        for i in 1:length(generated_ops)
            for j in i+1:length(generated_ops)
                if (i, j) ∉ checked_pair_indices
                    commutator = generated_ops[i] * generated_ops[j] - generated_ops[j] * generated_ops[i]
                    if real(tr(dagger(commutator) * commutator)) > tol
                        commutator = normalize_op(commutator)
                        if are_ops_linear_independent(vcat(generated_ops, [commutator]))
                            push!(generated_ops, commutator)
                        end
                        push!(checked_pair_indices, (i, j))
                    end
                end
            end
        end
    end
    return generated_ops
end

function is_controlable(ops_list)
    if length(ops_list) == 0
        return false
    end
    dim = length(ops_list[1].data)
    generated_ops = find_dynamical_lie_algebra(ops_list)
    return length(generated_ops) == dim - 1
end

is_controlable([Sx, Sz^2])

lie_algebra = find_dynamical_lie_algebra([Sz, Sx])

function find_adjoint_matrix(op, lie_algebra; tol=1e-10)
    n = length(lie_algebra)
    adjoint_matrix = zeros(ComplexF64, n, n)

    for i in 1:n
        commutator = op * lie_algebra[i] - lie_algebra[i] * op
        if real(tr(dagger(commutator) * commutator)) < tol
            adjoint_matrix[i, :] .= 0
        else
            for j in 1:n
                adjoint_matrix[i, j] = tr(dagger(lie_algebra[j]) * commutator)
            end
        end
    end

    return adjoint_matrix
end

function find_centralizer(op, lie_algebra; tol=1e-10)
    adjoint_matrix = find_adjoint_matrix(op, lie_algebra)
    U, s, V = svd(adjoint_matrix) # could use Krylov here instead?
    kernel_indices = findall(isapprox.(s, 0; atol=tol))
    centralizer_ops = [sum([V[i, j] * lie_algebra[j] for j in eachindex(lie_algebra)]) for i in kernel_indices]
    return centralizer_ops
end

find_adjoint_matrix(Sx, lie_algebra)

find_centralizer(Sx, lie_algebra)

function find_self_bracket(lie_algebra)
    generated_ops = []

    for i in eachindex(lie_algebra)
        for j in i+1:length(lie_algebra)
            commutator = lie_algebra[i] * lie_algebra[j] - lie_algebra[j] * lie_algebra[i]
            if real(tr(dagger(commutator) * commutator)) > tol
                commutator = normalize_op(commutator)
                if are_ops_linear_independent(vcat(lie_algebra, [commutator]))
                    push!(generated_ops, commutator)
                end
            end
        end
    end

    return generated_ops
end

function find_cartan_subalgebra(lie_algebra)

end

function render_dykin_diagram(lie_algebra)

end