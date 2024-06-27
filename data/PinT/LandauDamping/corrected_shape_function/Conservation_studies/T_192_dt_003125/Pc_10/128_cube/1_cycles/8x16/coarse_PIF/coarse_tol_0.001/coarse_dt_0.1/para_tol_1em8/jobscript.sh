#!/bin/bash -x
#SBATCH --account=cstma
#SBATCH --nodes=32
#SBATCH --ntasks-per-node=4
#SBATCH --output=landau.out
#SBATCH --error=landau.error
#SBATCH --time=00:10:00
#SBATCH --partition=booster
#SBATCH --gres=gpu:4
##SBATCH --gpus-per-task=1

#srun ./LandauDampingPinT 128 128 128 128 128 128 20971520 19.2 0.003125 0.1 1e-11 1 B-spline 1 8 16 0.001 1e-12 PIC --info 5
srun ./LandauDampingPinT 128 128 128 128 128 128 20971520 19.2 0.003125 0.1 1e-8 1 B-spline 1 8 16 0.001 1e-7 PIF --info 5
#srun ./LandauDampingPinT 128 128 128 128 128 128 20971520 19.2 0.05 0.1 1e-4 1 B-spline 1 8 16 0.001 1e-4 PIF --info 5