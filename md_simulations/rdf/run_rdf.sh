#!/bin/bash
#SBATCH --nodes=1
#SBATCH --partition=SGPU,GGPU
#SBATCH --gres=gpu:1 --ntasks=1
#SBATCH --time=24:00:00

source ~/.bashrc

conda activate omm8plumed

export OPENMM_CPU_THREADS=1
export OMP_NUM_THREADS=1

python run_rdf.py > rdf.out
