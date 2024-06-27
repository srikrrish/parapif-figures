#!/bin/bash -x
#SBATCH --account=cstma
#SBATCH --nodes=32
#SBATCH --ntasks-per-node=4
#SBATCH --output=penning.out
#SBATCH --error=penning.error
#SBATCH --time=00:10:00
#SBATCH --partition=booster
#SBATCH --gres=gpu:4
##SBATCH --gpus-per-task=1

#srun ./PenningTrapPinT 64 64 64 64 64 64 2621440 19.2 0.05 0.003125 1e-11 16 B-spline 1 4 32 0.1 1e-12 PIC --info 5
srun ./PenningTrapPinT 64 64 64 64 64 64 2621440 19.2 0.003125 0.003125 1e-8 16 B-spline 1 4 32 0.1 1e-7 PIC --info 5
#srun ./PenningTrapPinT 64 64 64 64 64 64 2621440 19.2 0.05 0.003125 1e-5 16 B-spline 1 4 32 0.1 1e-4 PIC --info 5