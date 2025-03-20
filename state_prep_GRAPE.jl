using QuantumOptics, Plots, Optim, FiniteDiff, Random, LaTeXStrings, Serialization, MLStyle, ArgParse
include("quantum_control/circuit_evolution.jl")
include("quantum_control/expm_utils.jl")

# Add Adam optimizer options to the argument parser
function my_parse_args()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--linesearch"
        help = "Line search method (e.g., BackTracking, HagerZhang, MoreThuente, Static)"
        arg_type = String
        default = ""
        "--results-dir"
        help = "Directory to save results"
        arg_type = String
        default = "results"
        "--g-tol"
        help = "Gradient tolerance"
        arg_type = Float64
        default = 1e-7
        "--m"
        help = "LBFGS memory parameter"
        arg_type = Int
        default = 20
        "--iterations"
        help = "Number of iterations"
        arg_type = Int
        default = 10
        "--alg-name"
        help = "Optimization algorithm name (e.g., LBFGS, GradientDescent, Adam)"
        arg_type = String
        default = "LBFGS"
        "--alpha"
        help = "Learning rate for Adam optimizer"
        arg_type = Float64
        default = 0.001
        "--beta-mean"
        help = "Exponential decay rate for the first moment estimates in Adam"
        arg_type = Float64
        default = 0.9
        "--beta-var"
        help = "Exponential decay rate for the second moment estimates in Adam"
        arg_type = Float64
        default = 0.999
        "--epsilon"
        help = "Small constant for numerical stability in Adam"
        arg_type = Float64
        default = 1e-8
    end
    return parse_args(s)
end

ARGS = my_parse_args()

for (arg, val) in ARGS
    println("  $arg  =>  $val")
end


RESULTS_DIR = ARGS["results-dir"]
alg_name = ARGS["alg-name"]
g_tol = ARGS["g-tol"]
m = ARGS["m"]
iterations = ARGS["iterations"]
alpha = ARGS["alpha"]
beta_mean = ARGS["beta-mean"]
beta_var = ARGS["beta-var"]
epsilon = ARGS["epsilon"]

# Define the cutoff for the Fock space
Nfock = 30

# Create the Fock space basis
b_fock = FockBasis(Nfock)
b_spin = SpinBasis(1 // 2)

# Create a composite Hilbert space
bases = [b_fock, b_spin]
full_basis = tensor(bases...)

# Function to tensor an operator with identity operators in all other bases
function mb(op, index, bases)
    ops = [index == i ? op : identityoperator(bases[i]) for i in eachindex(bases)]
    return tensor(ops...)
end

# Example usage of the mb function
a_op = mb(destroy(b_fock), 1, [b_fock, b_spin])
σ_minus = mb(sigmam(b_spin), 2, [b_fock, b_spin])


init_state = tensor(fockstate(b_fock, 0), spindown(b_spin))

targ_state = tensor((coherentstate(b_fock, 3.0) + coherentstate(b_fock, -3.0) + coherentstate(b_fock, 3.0im) + coherentstate(b_fock, -3.0im)), spindown(b_spin))
targ_state = targ_state / norm(targ_state)

linesearch = @match ARGS["linesearch"] begin
    "BackTracking" => Optim.LineSearches.BackTracking()
    "HagerZhang" => Optim.LineSearches.HagerZhang()
    "MoreThuente" => Optim.LineSearches.MoreThuente()
    "Static" => Optim.LineSearches.Static()
    "" => Optim.LineSearches.Static()
    _ => error("Unsupported linesearch type: $(ARGS["linesearch"])")
end

alg = @match alg_name begin
    "Adam" => Adam(; alpha=alpha,
        beta_mean=beta_mean,
        beta_var=beta_var,
        epsilon=epsilon)
    "AdaMax" => AdaMax(; alpha=alpha,
        beta_mean=beta_mean,
        beta_var=beta_var,
        epsilon=epsilon)
    "BFGS" => BFGS(; linesearch=linesearch)
    "LBFGS" => LBFGS(; m=m, linesearch=linesearch)

    "GradientDescent" => GradientDescent()
    _ => error("Unsupported algorithm: $alg_name")
end

U_q = Gate(params -> displacement(params, σ_minus / 2), params -> displacement_dparam(params, σ_minus / 2), nothing, 2, L"$U_q$")
U_qc = Gate(params -> displacement(params, dagger(a_op) * σ_minus / 2), params -> displacement_dparam(params, dagger(a_op) * σ_minus / 2), nothing, 2, L"$U_{qc}$")

targ_state = tensor((coherentstate(b_fock, 3.0) + coherentstate(b_fock, -3.0) + coherentstate(b_fock, 3.0im) + coherentstate(b_fock, -3.0im)), spindown(b_spin))
targ_state = targ_state / norm(targ_state)
fidelity_cat = Cost((params, ρ) -> 0, ρ -> 1 - real(tr(dagger(dm(targ_state)) * ρ)), (params, ρ) -> 0 * ρ, ρ -> -dagger(dm(targ_state)), (params, ρ) -> real.(0 * params))

UnitCircuit = QuantumCircuit(repeat([U_q, U_qc], outer=20))
params, result = GRAPE(dm(init_state), UnitCircuit, fidelity_cat; optim_alg=alg, optim_options=Optim.Options(g_tol=g_tol, show_trace=true, store_trace=true, iterations=iterations))

# Create the results directory if it doesn't exist
if !isdir(RESULTS_DIR)
    mkpath(RESULTS_DIR)
end

# Save the parameters, result, and command line arguments
params_file = joinpath(RESULTS_DIR, "params.jls")
result_file = joinpath(RESULTS_DIR, "result.jls")
args_file = joinpath(RESULTS_DIR, "args.jls")

serialize(params_file, params)
serialize(result_file, result)
serialize(args_file, ARGS)