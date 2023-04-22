#!/bin/bash
#SBATCH --partition=batch

#SBATCH --qos=240c-1h_batch
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=10G
#SBATCH --job-name="PC1"
#SBATCH --output=PC1.%J.out

#SBATCH --error=PC1.%J.err
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