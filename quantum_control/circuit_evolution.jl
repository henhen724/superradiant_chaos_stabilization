# module QuantumControl
using KernelAbstractions

import QuantumOptics: Operator

struct Gate
    func::Function #Vector{<:Number} -> Operator
    dfunc_dparam::Union{Function,Nothing} #Vector{<:Number} -> Vector{Operator}
    d2func_dparam2::Union{Function,Nothing} #Vector{<:Number} -> Matrix{Operator}
    num_params::Int
    name::Union{String,Nothing}
end

function Gate(H::Any; name=nothing)
    func = param -> exp(-im * param[1] * H)
    dfunc_dparam = param -> [-im * H * exp(-im * param[1] * H)]
    d2func_dparam2 = param -> [-H^2 * exp(-im * param[1] * H)]
    return Gate(func, dfunc_dparam, d2func_dparam2, 1, name)
end

struct Measurement
    func::Function #Vector{<:Number} -> Vector{Tuple{Operator,<:Number}}
    dfunc_dparam::Union{Function,Nothing} #Vector{<:Number} -> Vector{Vector{Tuple{Operator,<:Number}}}
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
    return timeseries, measurement_record
end

using Plots

function draw_circuit_element(element, x_position, param_str=nothing)
    if element isa Gate
        rect = Shape([x_position - 0.4, x_position - 0.4, x_position + 0.4, x_position + 0.4], [0.5, 1.5, 1.5, 0.5])
        plot!(rect, fillcolor=:blue, label="", seriestype=:shape, linejoin=:round)
        annotate!(x_position, 1.2, text(element.name, :center, 10, :white))
        if param_str !== nothing
            annotate!(x_position, 0.8, text("($param_str)", :center, 6, :white))
        end
    elseif element isa Measurement
        rect = Shape([x_position - 0.4, x_position - 0.4, x_position + 0.4, x_position + 0.4], [0.5, 1.5, 1.5, 0.5])
        plot!(rect, fillcolor=:green, label="", seriestype=:shape, linejoin=:round)
        annotate!(x_position, 1.2, text(element.name, :center, 10, :white))
        if param_str !== nothing
            annotate!(x_position, 0.8, text("($param_str)", :center, 6, :white))
        end
    else
        error("Unknown element type in quantum circuit")
    end
end

function draw_circuit(circuit::QuantumCircuit; params::Union{Nothing,Any}=nothing, fig=nothing)
    num_elements = length(circuit.elements)
    x_positions = 1:num_elements
    if fig === nothing
        fig = plot(legend=false, xlim=(0, num_elements + 1), ylim=(0, 2), aspect_ratio=:equal, framestyle=:none)
    else
        plot!(fig, legend=false, xlim=(0, num_elements + 1), ylim=(0, 2), aspect_ratio=:equal, framestyle=:none)
    end

    for (i, element) in enumerate(circuit.elements)
        if isa(params, Nothing)
            draw_circuit_element(element, x_positions[i])
        else
            param_str = join(round.(params[i], digits=1), ", ")
            draw_circuit_element(element, x_positions[i], param_str)
        end
        if i > 1
            plot!([x_positions[i-1] + 0.4, x_positions[i] - 0.4], [1, 1], color=:black, linewidth=2)
        end
    end

    # Draw time arrow
    arrow_x = [0.5, num_elements + 0.5]
    arrow_y = [0.2, 0.2]
    plot!(arrow_x, arrow_y, arrow=:arrow, color=:black, linewidth=2)
    annotate!((num_elements + 0.5) / 2, 0.1, text("Time", :center, 10, :black))

    return fig
end



struct Cost
    integrated_cost::Function # (param, ρ) -> <:Number
    terminal_cost::Function # ρ -> <:Number
    d_integrated_cost_dρ::Function # (param, ρ) -> Matrix
    d_terminal_cost_dρ::Function # ρ -> Matrix
    d_integrated_cost_dparams::Function # (param, ρ) -> Vector
end

using ForwardDiff, FiniteDiff

function Cost(integrated_cost::Function, terminal_cost::Function; autodiff=true)
    vec_integrated_cost = (param, ρ_vec, size) -> integrated_cost(param, reshape(ρ_vec, size))
    vec_terminal_cost = (ρ_vec, size) -> terminal_cost(reshape(ρ_vec, size))
    if autodiff
        d_integrated_cost_dρ = (param, ρ) -> reshape(ForwardDiff.gradient(ρ_vec -> vec_integrated_cost(param, ρ_vec, size(ρ)), reshape(dense(ρ.data), length(ρ.data))), size(ρ))
        d_terminal_cost_dρ = ρ -> reshape(ForwardDiff.gradient(ρ_vec -> vec_terminal_cost(param, ρ_vec, size(ρ)), reshape(dense(ρ.data), length(ρ.data))), size(ρ))
        d_integrated_cost_dparams = (param, ρ) -> ForwardDiff.gradient(param_vec -> integrated_cost(param_vec, ρ), param)
    end
    return Cost(integrated_cost, terminal_cost, d_integrated_cost_dρ, d_terminal_cost_dρ, d_integrated_cost_dparams)
end

function co_evolve(ρ_t, measurement_record, cost::Cost, circuit::QuantumCircuit, params::Vector{Vector{Float64}})
    σ_t = [0 * ρ for ρ in ρ_t]
    σ_t[end] = cost.d_terminal_cost_dρ(ρ_t[end])
    for (i, element) in reverse(collect(enumerate(circuit.elements)))
        param = params[i]
        if element isa Gate
            gate = element.func(param)
            σ_t[i] = dagger(gate) * σ_t[i+1] * gate
            # println("The costate norm", norm(σ_t[i]))
        elseif element isa Measurement
            Kraus_ops = element.func(params)
            println("TODO measurement costate not yet added")
            if averaged
                new_density_matrix = 0 * σ_t[i+1]
                for (op, label) in Kraus_ops
                    new_density_matrix += dagger(op) * σ_t[i+1] * op
                end
                σ_t[i] = new_density_matrix
            else
                selected_op, label = Kraus_ops[measurement_record[i][2]]
                σ_t[i] = dagger(selected_op) * σ_t[i+1] * selected_op / (probabilites[rand_index]) - dagger(selected_op) * selected_op * tr(dagger(selected_op) * σ_t[i+1] * selected_op * ρ_t[i]) / (probabilites[rand_index])^2
            end
        else
            error("Unknown element type in quantum circuit")
        end
        σ_t[i] = σ_t[i] + cost.d_integrated_cost_dρ(param, ρ_t[i+1])
    end
    return σ_t
end

# Kernel function to compute gradients for a single element
@kernel function compute_gradient_kernel(grad, @Const(ρ_t), @Const(σ_t), cost, circuit_elements, params)
    i = @index(Global)
    element = circuit_elements[i]
    param = params[i]
    if element isa Gate
        if isa(element.dfunc_dparam, Nothing)
            # Handle unimplemented case
            grad[i] .= 0.0 # Placeholder for unimplemented case
        else
            gate = element.func(param)
            dgate = element.dfunc_dparam(param)
            for j in eachindex(dgate)
                grad[i][j] = 2 * real(
                    tr(dagger(σ_t[i+1]) * gate * ρ_t[i] * dagger(dgate[j])) +
                    tr(dagger(σ_t[i+1]) * dgate[j] * ρ_t[i] * dagger(gate))
                )
            end
        end
    elseif element isa Measurement
        # Handle unimplemented Measurement case
        grad[i] .= 0.0 # Placeholder for unimplemented case
    else
        # Handle unknown element type
        grad[i] .= 0.0 # Placeholder for unknown element type
    end
end

function calculate_param_gradients(grad, ρ_t, σ_t, cost::Cost, circuit::QuantumCircuit, params::Vector{Vector{Float64}}; backend=CPU(), averaged=false)
    # Prepare data for the kernel
    circuit_elements = circuit.elements
    num_elements = length(circuit_elements)

    # Launch the kernel
    compute_gradient_kernel(backend, 20)(grad, ρ_t, σ_t, cost, circuit_elements, params, ndrange=size(circuit_elements))
    synchronize(backend)
end

function reshape_to_params_vector(flat_params, num_params_vec)
    params::Vector{Vector{eltype(flat_params)}} = []
    num_of_params_so_far = 0
    for num_params in num_params_vec
        push!(params, view(flat_params, num_of_params_so_far+1:num_of_params_so_far+num_params))
        num_of_params_so_far += num_params
    end
    return params
end

function clip_val(val, min_v, max_v)
    return min(max(val, min_v), max_v)
end

function GRAPE(ρ_0, circuit::QuantumCircuit, cost::Cost; seed=nothing, optim_alg=LBFGS(; m=20), optim_options=Optim.Options(g_tol=1e-7, show_trace=true, store_trace=true, iterations=1000))
    Random.seed!(seed)
    inital_params = [2 * π * randn(Float64, ele.num_params) for ele in circuit.elements]
    num_params_vec = [ele.num_params for ele in circuit.elements]
    function objective(x)
        params = reshape_to_params_vector(x, num_params_vec)
        ρ_t, _ = evolve(ρ_0, circuit, params)
        integrated_cost_sum = sum(cost.integrated_cost(params[i], ρ) for (i, ρ) in enumerate(ρ_t[begin:end-1]))
        terminal_cost_value = cost.terminal_cost(ρ_t[end])
        return integrated_cost_sum + terminal_cost_value
    end

    function gradient!(G, x)
        params = reshape_to_params_vector(x, num_params_vec)
        ρ_t, _ = measurement_record = evolve(ρ_0, circuit, params)
        σ_t = co_evolve(ρ_t, measurement_record, cost, circuit, params)
        grad = [zeros(Float64, ele.num_params) for ele in circuit.elements]
        calculate_param_gradients(grad, ρ_t, σ_t, cost, circuit, params; backend=CPU())
        G .= vcat(grad...)#map(x -> clip_val(x, -10.0, 10.0), vcat(grad...))
    end

    x0 = vcat(inital_params...)
    result = optimize(objective, gradient!, x0, optim_alg, optim_options)

    return reshape_to_params_vector(result.minimizer, num_params_vec), result
end

# end # module QuantumControl