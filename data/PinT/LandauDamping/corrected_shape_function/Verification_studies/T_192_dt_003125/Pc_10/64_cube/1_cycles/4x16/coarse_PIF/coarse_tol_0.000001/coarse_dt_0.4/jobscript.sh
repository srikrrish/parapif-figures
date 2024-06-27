#!/bin/bash -x
#SBATCH --account=cstma
#SBATCH --nodes=16
#SBATCH --ntasks-per-node=4
#SBATCH --output=landau.out
#SBATCH --error=landau.error
#SBATCH --time=00:15:00
#SBATCH --partition=booster
#SBATCH --gres=gpu:4
##SBATCH --gpus-per-task=1

srun ./LandauDampingPinT 64 64 64 64 64 64 2621440 19.2 0.003125 0.4 1e-11 1 B-spline 1 4 16 0.000001 1e-6 PIF --info 5
#srun ./LandauDampingPinT 64 64 64 64 64 64 2621440 19.2 0.05 0.4 1e-5 1 B-spline 1 4 16 0.000001 1e-4 PIC --info 5
