#!/bin/bash
#SBATCH --partition=serial

#SBATCH --qos=84c-1d_serial
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=10G
#SBATCH --job-name="MPC9_1"
#SBATCH --output=MPC9_1.%J.out

#SBATCH --error=MPC9_1.%J.err
#SBATCH --mail-type=ALL
#SBATCH --requeue

echo "SLURM_JOBID="$SLURM_JOBID
echo "SLURM_JOB_NODELIST"=$SLURM_JOB_NODELIST
echo "SLURM_NNODES"=$SLURM_NNODES
echo "SLURMTMPDIR="$SLURMTMPDIR
echo "working directory = "$SLURM_SUBMIT_DIR

ulimit -s unlimited

# MAIN
bash commands.sh