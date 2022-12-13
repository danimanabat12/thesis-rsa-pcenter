#! /bin/bash
for d in */*/; do
    # Will print */ if no directories are available
    STR=$(pwd)
    cd $(pwd)/$d
    sbatch compile.slurm
    cd $STR
    # sbatch $d
done
