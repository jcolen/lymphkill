#!/bin/bash
#SBATCH --time=1-00:00:00
#SBATCH --output=mask_generation_proton.out
#SBATCH --error=mask_generation_proton.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jc8kf@virginia.edu
#SBATCH --account=phys_nrf
#SBATCH --partition=standard
#SBATCH --mem=30000
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --tasks=1

#./mask_generation_proton.sh BM BM_organs_short.mat clinical/CT_1.2.840.113704.1.111.588.1370268691.469.dcm BM_proton_RP_RD/RD1.2.752.243.1.1.20180911152103617.3700.62580.dcm BM_masks_proton.mat
#./mask_generation_proton.sh ML ML_organs_new.mat ML_Clinical/CT_1.2.840.113704.1.111.3676.1409664164.575.dcm ML_proton_RP_RD/RD1.2.752.243.1.1.20180920171152288.2200.13045.dcm ML_masks_proton.mat
./mask_generation_proton.sh LL LL_organs_short.mat LL/CT_1.2.840.113704.1.111.3520.1441734056.917.dcm LL_proton/RD1.2.752.243.1.1.20180927103209566.1000.60835.dcm LL_masks_proton.mat

