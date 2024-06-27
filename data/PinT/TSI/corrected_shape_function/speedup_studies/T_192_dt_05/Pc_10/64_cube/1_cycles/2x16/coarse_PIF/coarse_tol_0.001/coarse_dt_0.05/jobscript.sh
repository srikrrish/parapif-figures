#!/bin/bash -x
#SBATCH --account=cstma
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=4
#SBATCH --output=tsi.out
#SBATCH --error=tsi.error
#SBATCH --time=00:10:00
#SBATCH --partition=booster
#SBATCH --gres=gpu:4
##SBATCH --gpus-per-task=1

#srun ./BumponTailInstabilityPinT 64 64 64 64 64 64 2621440 19.2 0.003125 0.05 1e-8 1 B-spline 1 2 16 0.001 1e-7 PIF --info 5
srun ./BumponTailInstabilityPinT 64 64 64 64 64 64 2621440 19.2 0.05 0.05 1e-5 1 B-spline 1 2 16 0.001 1e-4 PIF --info 5
