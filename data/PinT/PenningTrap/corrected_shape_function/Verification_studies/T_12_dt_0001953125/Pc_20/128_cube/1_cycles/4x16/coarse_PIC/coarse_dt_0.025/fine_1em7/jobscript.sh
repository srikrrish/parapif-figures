#!/bin/bash -x
#SBATCH --account=cstma
#SBATCH --nodes=16
#SBATCH --ntasks-per-node=4
#SBATCH --output=penning.out
#SBATCH --error=penning.error
#SBATCH --time=00:15:00
#SBATCH --partition=booster
#SBATCH --gres=gpu:4
##SBATCH --gpus-per-task=1

srun ./PenningTrapPinT 128 128 128 128 128 128 41943040 1.2 0.0001953125 0.025 1e-11 1 B-spline 1 4 16 0.1 1e-7 PIC --info 5
#srun ./PenningTrapPinT 128 128 128 128 128 128 41943040 19.2 0.003125 0.025 1e-8 1 B-spline 7 4 16 0.1 1e-7 PIC --info 5
#srun ./PenningTrapPinT 128 128 128 128 128 128 41943040 19.2 0.05 0.025 1e-5 1 B-spline 1 4 16 0.1 1e-4 PIC --info 5
