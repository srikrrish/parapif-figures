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

#srun ./PenningTrapPinT 64 64 64 64 64 64 5242880 1.2 0.0001953125 0.00009765625 1e-11 1 B-spline 1 4 16 0.000001 1e-6 PIF --info 5
srun ./PenningTrapPinT 64 64 64 64 64 64 5242880 0.075 0.00001220703125 0.00009765625 1e-11 1 B-spline 1 4 16 0.000001 1e-6 PIF --info 5
#srun ./PenningTrapPinT 64 64 64 64 64 64 5242880 19.2 0.003125 0.00009765625 1e-8 1 B-spline 7 4 16 0.000001 1e-7 PIC --info 5
#srun ./PenningTrapPinT 64 64 64 64 64 64 5242880 19.2 0.05 0.00009765625 1e-5 1 B-spline 1 4 16 0.000001 1e-4 PIC --info 5
