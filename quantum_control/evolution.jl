# module QuantumControl

import QuantumOptics: Operator

struct Gate
    func::Function #Vector{<:Number} -> Vector{Tuple{Operator,<:Number}}
    num_params::Int
    name::Union{String,Nothing}
end


struct Measurement
    func::Function #Vector{<:Number} -> Vector{Tuple{Operator,<:Number}}
    num_params::Int
    name::Union{String,Nothing}
end

struct QuantumCircuit
    elements::Vector{Union{Gate,Measurement}}
end

# Define the evolve function
function evolve(density_matrix, circuit::QuantumCircuit, params::Vector{Vector{Float64}}; averaged=false, record_timeseries=true)
    timeseries = []
    measurement_record = []
    if record_timeseries
        push!(timeseries, copy(density_matrix))
    end
    for (i, element) in enumerate(circuit.elements)
        param = params[i]
        if element isa Gate
            gate = element.func(param)
            density_matrix = gate * density_matrix * dagger(gate)
        elseif element isa Measurement
            Kraus_ops = element.func(params)
            if averaged
                new_density_matrix = 0 * density_matrix
                for (op, label) in Kraus_ops
                    new_density_matrix += op * density_matrix * dagger(op)
                end
                density_matrix = new_density_matrix
            else
                probabilites = [tr(op * density_matrix * dagger(op)) for (op, label) in Kraus_ops]
                uni_rand = rand()
                rand_index = findfirst(x -> x >= uni_rand, cumsum(probabilites))
                selected_op, label = Kraus_ops[rand_index]
                density_matrix = selected_op * density_matrix * dagger(selected_op) / (probabilites[rand_index])
                push!(measurement_record, (element.name, label, rand_index))
            end
        else
            error("Unknown element type in quantum circuit")
        end
        if record_timeseries
            push!(timeseries, copy(density_matrix))
        end
    end
    return density_matrix, timeseries, measurement_record
end

struct Cost
    integrated_cost::Function
    terminal_cost::Function
    grad_integrated_cost::Function
    grad_terminal_cost::Function
end

function reshape_to_params_vector(flat_params, num_params_vec)
    params = []
    num_of_params_so_far = 0
    for num_params in num_params_vec
        push!(params, view(flat_params, num_of_params_so_far:num_of_params_so_far+num_params))
        num_of_params_so_far += num_params
    end
    return params
end

function GRAPE(density_matrix, circuit::QuantumCircuit, cost::Cost; seed=nothing)
    Random.seed!(seed)
    inital_params = [2 * π * randn(ele.num_params) for ele in circuit.elements]
    function objective(x)
        params = reshape_to_params_vector(x, [ele.num_params for ele in circuit.elements])
        ρ_t, _, _ = evolve(density_matrix, circuit, params)
        integrated_cost_sum = sum(cost.integrated_cost(ρ) for ρ in ρ_t)
        terminal_cost_value = cost.terminal_cost(ρ_t[end])
        return integrated_cost_sum + terminal_cost_value
    end
    end

    function gradient!(G, x)
        ρ_t, _, measurement_record = evolve(density_matrix, circuit, reshape_to_params_vector(x, [ele.num_params for ele in circuit.elements]))
        σ_t = co_evolve(ρ_t, measurement_record, cost, circuit, reshape_to_params_vector(x, [ele.num_params for ele in circuit.elements]))
        grad = reshape_to_params_vector(G, [ele.num_params for ele in circuit.elements])
        for (i, element) in enumerate(circuit.elements)
            param = reshape_to_params_vector(x, [ele.num_params for ele in circuit.elements])[i]
            if element isa Gate
                gate = element.func(param)
                grad[i] .= 2 * real.(tr(dagger(gate) * σ_t[i+1] * gate * ρ_t[i] - σ_t[i+1] * gate * ρ_t[i] * dagger(gate)))
            elseif element isa Measurement
                Kraus_ops = element.func(param)
                if averaged
                    for (op, label) in Kraus_ops
                        grad[i] .= grad[i] .+ 2 * real.(tr(dagger(op) * σ_t[i+1] * op * ρ_t[i] - σ_t[i+1] * op * ρ_t[i] * dagger(op)))
                    end
                else
                    selected_op, label = Kraus_ops[measurement_record[i][2]]
                    grad[i] .= 2 * real.(tr(dagger(selected_op) * σ_t[i+1] * selected_op * ρ_t[i] - σ_t[i+1] * selected_op * ρ_t[i] * dagger(selected_op)))
                end
            else
                error("Unknown element type in quantum circuit")
            end
        end
    end

    x0 = vcat(inital_params...)
    result = optimize(objective, gradient!, x0, LBFGS(), Optim.Options(g_tol=1e-5, show_trace=true, iterations=1000))

    result.minimizer
end

function co_evolve(ρ_t, measurement_record, cost::Cost, circuit::QuantumCircuit, params::Vector{Vector{Float64}})
    σ_t = [0 * ρ for ρ in ρ_t]
    σ_t[end] = cost.grad_terminal_cost(ρ_t[end])
    for (i, element) in reverse(collect(enumerate(circuit.elements)))
        param = params[i]
        if element isa Gate
            gate = element.func(param)
            σ_t[i] = dagger(gate) * σ_t[i+1] * gate + cost.grad_integrated_cost(ρ_t[i+1])
        elseif element isa Measurement
            Kraus_ops = element.func(params)
            println("TODO measurement costate not yet added")
            if averaged
                new_density_matrix = 0 * σ_t[i+1]
                for (op, label) in Kraus_ops
                    new_density_matrix += dagger(op) * σ_t[i+1] * op
                end
                σ_t[i] = new_density_matrix + cost.grad_integrated_cost(ρ_t[i+1])
            else
                selected_op, label = findfirst(Kraus_ops[measurement_record[i][2]])
                σ_t[i] = dagger(selected_op) * σ_t[i+1] * selected_op / (probabilites[rand_index]) - dagger(selected_op) * selected_op * tr(dagger(selected_op) * σ_t[i+1] * selected_op * ρ_t[i]) / (probabilites[rand_index])^2
            end
        else
            error("Unknown element type in quantum circuit")
        end
    end
    return σ_t
end

# end # module QuantumControl