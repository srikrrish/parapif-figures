#!/bin/bash -x
#SBATCH --account=cstma
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --output=tsi.out
#SBATCH --error=tsi.err
#SBATCH --time=00:10:00
#SBATCH --partition=booster
#SBATCH --gres=gpu:4
##SBATCH --gpus-per-task=1

#srun ./LandauDamping 512 512 512 1073741824 20 FFT 1.0 2.0 --info 10

#/p/software/juwelsbooster/stages/2023/software/Nsight-Systems/2023.2.1-GCCcore-11.3.0/bin/nsys profile --trace=mpi,cuda,nvtx --stats=true --force-overwrite=true srun ./LandauDampingPIF 64 64 64 2621440 6144 0.05 B-spline 1 --info 5
srun ./BumponTailInstabilityPIF 64 64 64 2621440 6144 0.003125 B-spline 1 1e-7 --info 5
