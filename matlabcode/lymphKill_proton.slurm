#!/bin/bash
#SBATCH --time=5-00:00:00
#SBATCH --output=lymphKill_proton.out
#SBATCH --error=lymphKill_proton.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jc8kf@virginia.edu
#SBATCH --account=phys_nrf
#SBATCH --partition=standard
#SBATCH --mem=30000
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --tasks=1

#./lymphKill_proton.sh BM BM_masks_proton.mat 60.0 0 4 BM_proton_RP_RD_filenames.txt BM/BM_proton_RP_RD_dose.mat
#./lymphKill_proton.sh ML ML_masks_proton.mat 60.0 0 4 ML_proton_RP_RD_filenames.txt ML/ML_proton_RP_RD_dose.mat
./lymphKill_proton.sh LL LL_masks_proton.mat 60.0 0 4 LL_proton_filenames.txt LL/LL_proton_dose.mat

