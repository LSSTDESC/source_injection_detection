#!/bin/bash
#SBATCH -A m1727
#SBATCH -C gpu
#SBATCH -q regular
#SBATCH -t 02:00:00
#SBATCH -c 8
#SBATCH -J gen1y
#SBATCH -o logs/slsim_%A_%a.out
#SBATCH -e logs/slsim_%A_%a.err

# -------------------------
# to call this: sbatch generate_slsim_auto.sh 1000 3 1 3000sqdeg_lsst_1y_sample_1.h5
# where you sample over 1000 sq deg 3 times, and create 1 year baseline light curves and save to 3000sqdeg_lsst_1y_sample_1.h5 
# -------------------------
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

module load conda
conda activate /global/common/software/m1727/vpadma/slsim_env

python generate_slsim_lsst_dataset.py --sky-area $1 --n-iter $2 --baseline $3 --output-file /pscratch/sd/v/vpadma/lens_finding/$4