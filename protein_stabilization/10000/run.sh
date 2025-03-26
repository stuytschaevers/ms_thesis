#!/bin/bash -l
#
#SBATCH --job-name=protStablization
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=1000
#SBATCH --partition=tier3
#SBATCH --gres=gpu:a100:1
#SBATCH --account=coin

echo "Current Folder:"
pwd
echo $PWD

spack load /l2gshlg

echo "Loaded Modules"
spack find --loaded

echo "Running simulation:"
python stab_rc_new.py
