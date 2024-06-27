#!/bin/bash -x
#SBATCH --account=cstma
#SBATCH --nodes=64
#SBATCH --ntasks-per-node=4
#SBATCH --output=penning.out
#SBATCH --error=penning.error
#SBATCH --time=00:15:00
#SBATCH --partition=booster
#SBATCH --gres=gpu:4
##SBATCH --gpus-per-task=1

srun ./PenningTrapPinT 16 16 16 16 16 16 40960 19.2 0.05 0.05 1e-11 1 B-spline 1 16 16 0.000001 1e-12 PIC --info 5
#srun ./PenningTrapPinT 16 16 16 16 16 16 40960 19.2 0.003125 0.05 1e-9 1 B-spline 1 16 16 0.000001 1e-7 PIC --info 5
#srun ./PenningTrapPinT 16 16 16 16 16 16 40960 19.2 0.05 0.05 1e-5 1 B-spline 1 16 16 0.000001 1e-4 PIC --info 5
