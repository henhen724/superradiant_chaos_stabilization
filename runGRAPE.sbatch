#!/bin/bash
#SBATCH --job-name=GRAPE
#SBATCH --output=grape_%A_%a.out
#SBATCH --error=grape_%A_%a.err
#SBATCH --array=0-9
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --mem=4G

# Load Julia module
module load julia

# Define the root directory for results
RESULTS_ROOT="results"

# Define the parameter combinations (no more than 10)
declare -a PARAM_COMBINATIONS=(
    "LBFGS,,1e-7,20,100"
    "LBFGS,BackTracking,1e-6,15,200"
    "Adam,,,0.001,0.9,0.999,1e-8"
    "Adam,,,0.01,0.8,0.99,1e-7"
    "AdaMax,,,0.001,0.9,0.999,1e-8"
    "AdaMax,,,0.01,0.8,0.99,1e-7"
    "GradientDescent,,,0.01"
    "GradientDescent,,,0.001"
    "BFGS,MoreThuente,1e-7"
    "BFGS,HagerZhang,1e-6"
)

# Get the parameter combination for this array task
PARAMS=${PARAM_COMBINATIONS[$SLURM_ARRAY_TASK_ID]}

# Parse the parameters
IFS=',' read -r alg_name linesearch g_tol m iterations alpha beta_mean beta_var epsilon <<< "$PARAMS"

# Run 10 iterations for each parameter combination
for run_id in {1..10}; do
    # Create a unique results directory for this run
    RESULTS_DIR="${RESULTS_ROOT}/${alg_name}_${linesearch}_${g_tol}_${m}_${iterations}_${alpha}_${beta_mean}_${beta_var}_${epsilon}/run_${run_id}"
    mkdir -p "$RESULTS_DIR"

    # Run the Julia script with the specified parameters
    julia state_prep_GRAPE.jl \
        --alg-name "$alg_name" \
        --linesearch "$linesearch" \
        --g-tol "$g_tol" \
        --m "$m" \
        --iterations "$iterations" \
        --alpha "$alpha" \
        --beta-mean "$beta_mean" \
        --beta-var "$beta_var" \
        --epsilon "$epsilon" \
        --results-dir "$RESULTS_DIR"
done