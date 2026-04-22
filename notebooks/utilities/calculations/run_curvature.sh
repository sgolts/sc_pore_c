#!/bin/bash
#SBATCH --job-name=curvature
#SBATCH --account=indikar1
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=logs/curvature_%j.out
#SBATCH --error=logs/curvature_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=sgolts@umich.edu

mkdir -p logs

CONDA_BASE=$(conda info --base)
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate /scratch/indikar_root/indikar1/sgolts/conda-envs/pore_c

INPUT="/nfs/turbo/umms-indikar/shared/projects/poreC/pipeline_outputs/higher_order/anndata/population_mESC_1000000_features.h5ad"
EDGE_OUT="outputs/edges_curvature.csv"
NODE_OUT="outputs/nodes_curvature.csv"

mkdir -p outputs

/scratch/indikar_root/indikar1/sgolts/conda-envs/pore_c/bin/python curvature.py \
    "$INPUT" \
    "$EDGE_OUT" \
    "$NODE_OUT" \