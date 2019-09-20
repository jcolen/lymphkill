#!/bin/bash
#SBATCH --time=5-00:00:00
#SBATCH --output=printbloodVelocities.out
#SBATCH --error=printbloodVelocities.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jc8kf@virginia.edu
#SBATCH --account=phys_nrf
#SBATCH --partition=standard
#SBATCH --mem=30000
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --tasks=1

./printBloodVelocities.sh BD BD_masks_fix.mat 29.3 0 5 BD_clinical_filenames.txt fix/BD_clinical_fix.mat
./printBloodVelocities.sh BF BF_masks_fix.mat 29.0 0 5 BF_clinical_filenames.tx fix/BF_clinical_fix.mat
./printBloodVelocities.sh SJ-1A SJ-1A_masks_fix.mat 33.3 1 5 SJ-1A_clinical_filenames.tx fix/SJ-1A_clinical_fix.mat
./printBloodVelocities.sh JL1 JL1_masks_fix.mat 50.9 0 5 JL1_clinical_filenames.txt fix/JL1_clinical_fix.mat
./printBloodVelocities.sh SR_1A SR_1A_masks_fix.mat 26.7 0 5 SR_clinical_filenames.txt fix/SR_1A_clinical_fix.mat
./printBloodVelocities.sh BB BB_masks_fix.mat 76.4 0 5 BB_clinical_filenames.txt fix/BB_clinical_fix.mat
./printBloodVelocities.sh BC BC_masks_fix.mat 20.0 0 5 BC_clinical_filenames.txt fix/BC_clinical_fix.mat
./printBloodVelocities.sh HR HR_masks_fix.mat 22.7 1 5 HR_clinical_filenames.txt fix/HR_clinical_fix.mat
./printBloodVelocities.sh BM BM_masks_fix.mat 24.7 1 5 BM_clinical_filenames.txt fix/BM_clinical_fix.mat
./printBloodVelocities.sh MM1 MM1_masks_fix.mat 31.7 0 5 MM1_clinical_filenames.txt fix/MM1_clinical_fix.mat
./printBloodVelocities.sh WJ WJ_masks_fix.mat 47.9 0 5 WJ_clinical_filenames.txt fix/WJ_clinical_fix.mat
./printBloodVelocities.sh BJ BJ_masks_fix.mat 21.9 0 5 BJ_clinical_filenames.txt fix/BJ_clinical_fix.mat
./printBloodVelocities.sh WH WH_masks_fix.mat 25.6 0 5 WH_clinical_filenames.txt fix/WH_clinical_fix.mat
./printBloodVelocities.sh JL JL_masks_fix.mat 22.5 1 5 JL_clinical_filenames.txt fix/JL_clinical_fix.mat
./printBloodVelocities.sh KE-1A KE-1A_masks_fix.mat 42.6 0 5 KE-1A_clinical_filenames.txt fix/KE-1A_clinical_fix.mat
./printBloodVelocities.sh HI HI_masks_fix.mat 136.9 0 5 HI_clinical_filenames.txt fix/HI_clinical_fix.mat
./printBloodVelocities.sh ME1 ME1_masks_fix.mat 29.2 1 5 ME1_clinical_filenames.tx fix/ME1_clinical_fix.mat
./printBloodVelocities.sh HH-1BResim HH-1BResim_masks_fix.mat 43.4 1 5 HH-1BResim_clinical_filenames.tx fix/HH-1BResim_clinical_fix.mat
./printBloodVelocities.sh BR BR_masks_fix.mat 24.9 0 5 BR_clinical_filenames.txt fix/BR_clinical_fix.mat
./printBloodVelocities.sh LV-1A LV-1A_masks_fix.mat 34.8 0 5 LV-1A_clinical_filenames.tx fix/LV-1A_clinical_fix.mat
./printBloodVelocities.sh HB HB_masks_fix.mat 35.4 0 5 HB_clinical_filenames.tx fix/HB_clinical_fix.mat
./printBloodVelocities.sh WC WC_masks_fix.mat 20.0 0 5 WC_clinical_filenames.txt fix/WC_clinical_fix.mat
./printBloodVelocities.sh FJ FJ_masks_fix.mat 16.9 0 5 FJ_clinical_filenames.txt fix/FJ_clinical_fix.mat
./printBloodVelocities.sh HA HA_masks_fix.mat 74.0 0 5 HA_clinical_filenames.txt fix/HA_clinical_fix.mat
./printBloodVelocities.sh JF JF_masks_fix.mat 28.1 0 5 JF_clinical_filenames.txt fix/JF_clinical_fix.mat
./printBloodVelocities.sh HD HD_masks_fix.mat 24.3 0 5 HD_clinical_filenames.tx fix/HD_clinical_fix.mat
./printBloodVelocities.sh LD LD_masks_fix.mat 20.7 1 5 LD_clinical_filenames.tx fix/LD_clinical_fix.mat
./printBloodVelocities.sh SJ2 SJ2_masks_fix.mat 55.2 0 5 SJ2_clinical_filenames.txt fix/SJ2_clinical_fix.mat
./printBloodVelocities.sh SR1 SR1_masks_fix.mat 53.6 0 5 SR1_clinical_filenames.txt fix/SR1_clinical_fix.mat
./printBloodVelocities.sh LL LL_masks_fix.mat 41.2 0 5 LL_clinical_filenames.txt fix/LL_clinical_fix.mat
./printBloodVelocities.sh SB-1A SB-1A_masks_fix.mat 27.7 0 5 SB-1A_clinical_filenames.tx fix/SB-1A_clinical_fix.mat
./printBloodVelocities.sh HR1 HR1_masks_fix.mat 32.0 0 5 HR1_clinical_filenames.txt fix/HR1_clinical_fix.mat
./printBloodVelocities.sh BL BL_masks_fix.mat 62.5 0 5 BL_clinical_filenames.txt fix/BL_clinical_fix.mat
./printBloodVelocities.sh RD-1A RD-1A_masks_fix.mat 37.7 0 5 RD-1A_clinical_filenames.txt fix/RD-1A_clinical_fix.mat
./printBloodVelocities.sh SH SH_masks_fix.mat 18.5 0 5 SH_clinical_filenames.txt fix/SH_clinical_fix.mat
./printBloodVelocities.sh BD1 BD1_masks_fix.mat 38.4 0 5 BD1_clinical_filenames.txt fix/BD1_clinical_fix.mat
./printBloodVelocities.sh MJ1 MJ1_masks_fix.mat 17.9 1 5 MJ1_clinical_filenames.txt fix/MJ1_clinical_fix.mat
./printBloodVelocities.sh PJ PJ_masks_fix.mat 14.3 0 5 PJ_clinical_filenames.txt fix/PJ_clinical_fix.mat
./printBloodVelocities.sh RM RM_masks_fix.mat 34.3 0 5 RM_clinical_filenames.txt fix/RM_clinical_fix.mat
./printBloodVelocities.sh PH PH_masks_fix.mat 61.1 0 5 PH_clinical_filenames.txt fix/PH_clinical_fix.mat
./printBloodVelocities.sh RH RH_masks_fix.mat 56.1 0 5 RH_clinical_filenames.tx fix/RH_clinical_fix.mat
./printBloodVelocities.sh VJ VJ_masks_fix.mat 23.2 0 5 VJ_clinical_filenames.txt fix/VJ_clinical_fix.mat
./printBloodVelocities.sh BB1 BB1_masks_fix.mat 51.1 0 5 BB1_clinical_filenames.txt fix/BB1_clinical_fix.mat
./printBloodVelocities.sh RA RA_masks_fix.mat 21.5 0 5 RA_clinical_filenames.txt fix/RA_clinical_fix.mat
./printBloodVelocities.sh AA AA_masks_fix.mat 46.3 0 5 AA_clinical_filenames.txt fix/AA_clinical_fix.mat
./printBloodVelocities.sh GV GV_masks_fix.mat 26.7 1 5 GV_clinical_filenames.txt fix/GV_clinical_fix.mat
./printBloodVelocities.sh AP AP_masks_fix.mat 28.4 0 5 AP_clinical_filenames.tx fix/AP_clinical_fix.mat
./printBloodVelocities.sh GS GS_masks_fix.mat 26.9 1 5 GS_clinical_filenames.txt fix/GS_clinical_fix.mat
./printBloodVelocities.sh CM CM_masks_fix.mat 32.6 1 5 CM_clinical_filenames.txt fix/CM_clinical_fix.mat
./printBloodVelocities.sh CL_1A CL_1A_masks_fix.mat 25.2 0 5 CL_1A_clinical_filenames.txt fix/CL_1A_clinical_fix.mat
./printBloodVelocities.sh JS-1B JS-1B_masks_fix.mat 22.6 0 5 JS-1B_clinical_filenames.txt fix/JS-1B_clinical_fix.mat
./printBloodVelocities.sh FL-1A FL-1A_masks_fix.mat 52.2 0 5 FL-1A_clinical_filenames.tx fix/FL-1A_clinical_fix.mat
./printBloodVelocities.sh CE CE_masks_fix.mat 38.6 0 5 CE_clinical_filenames.txt fix/CE_clinical_fix.mat
./printBloodVelocities.sh MJ MJ_masks_fix.mat 25.2 0 5 MJ_clinical_filenames.txt fix/MJ_clinical_fix.mat
./printBloodVelocities.sh GF GF_masks_fix.mat 35.0 0 5 GF_clinical_filenames.txt fix/GF_clinical_fix.mat
./printBloodVelocities.sh GE GE_masks_fix.mat 40.1 0 5 GE_clinical_filenames.txt fix/GE_clinical_fix.mat
./printBloodVelocities.sh AR AR_masks_fix.mat 30.1 0 5 AR_clinical_filenames.txt fix/AR_clinical_fix.mat
./printBloodVelocities.sh EC EC_masks_fix.mat 64.9 0 5 EC_clinical_filenames.txt fix/EC_clinical_fix.mat
./printBloodVelocities.sh LR_1A LR_1A_masks_fix.mat 28.8 0 5 LR_1A_clinical_filenames.txt fix/LR_1A_clinical_fix.mat
./printBloodVelocities.sh CS CS_masks_fix.mat 66.9 0 5 CS_clinical_filenames.tx fix/CS_clinical_fix.mat
./printBloodVelocities.sh GJ GJ_masks_fix.mat 35.4 0 5 GJ_clinical_filenames.tx fix/GJ_clinical_fix.mat
./printBloodVelocities.sh WJ1 WJ1_masks_fix.mat 40.3 0 5 WJ1_clinical_filenames.txt fix/WJ1_clinical_fix.mat
./printBloodVelocities.sh ME ME_masks_fix.mat 38.8 0 5 ME_clinical_filenames.txt fix/ME_clinical_fix.mat
./printBloodVelocities.sh HJ HJ_masks_fix.mat 53.2 0 5 HJ_clinical_filenames.txt fix/HJ_clinical_fix.mat
./printBloodVelocities.sh JM-1A JM-1A_masks_fix.mat 38.3 1 3 JM-1A_clinical_filenames.txt fix/JM-1A_clinical_fix.mat
./printBloodVelocities.sh ST ST_masks_fix.mat 32.6 1 5 ST_clinical_filenames.tx fix/ST_clinical_fix.mat
./printBloodVelocities.sh OH OH_masks_fix.mat 49.3 0 5 OH_clinical_filenames.txt fix/OH_clinical_fix.mat
./printBloodVelocities.sh MM MM_masks_fix.mat 49.3 0 5 MM_clinical_filenames.txt fix/MM_clinical_fix.mat
./printBloodVelocities.sh CA-2ALUL CA-2ALUL_masks_fix.mat 58.3 0 5 CA-2ALUL_clinical_filenames.txt fix/CA-2ALUL_clinical_fix.mat
./printBloodVelocities.sh MB MB_masks_fix.mat 40.4 0 3 MB_clinical_filenames.tx fix/MB_clinical_fix.mat
./printBloodVelocities.sh HM-1A HM-1A_masks_fix.mat 29.3 1 5 HM-1A_clinical_filenames.tx fix/HM-1A_clinical_fix.mat
./printBloodVelocities.sh SJ1-1A SJ1-1A_masks_fix.mat 40.5 1 5 SJ1-1A_clinical_filenames.txt fix/SJ1-1A_clinical_fix.mat
./printBloodVelocities.sh MK MK_masks_fix.mat 18.6 0 5 MK_clinical_filenames.txt fix/MK_clinical_fix.mat
./printBloodVelocities.sh TC TC_masks_fix.mat 26.2 0 5 TC_clinical_filenames.txt fix/TC_clinical_fix.mat
./printBloodVelocities.sh CE1 CE1_masks_fix.mat 28.2 0 5 CE1_clinical_filenames.txt fix/CE1_clinical_fix.mat
./printBloodVelocities.sh MD-1A MD-1A_masks_fix.mat 31.9 0 3 MD-1A_clinical_filenames.txt fix/MD-1A_clinical_fix.mat
./printBloodVelocities.sh MS MS_masks_fix.mat 47.1 0 5 MS_clinical_filenames.txt fix/MS_clinical_fix.mat
./printBloodVelocities.sh MR MR_masks_fix.mat 42.5 0 5 MR_clinical_filenames.tx fix/MR_clinical_fix.mat
./printBloodVelocities.sh CJ CJ_masks_fix.mat 31.4 1 5 CJ_clinical_filenames.tx fix/CJ_clinical_fix.mat
./printBloodVelocities.sh OP-1A OP-1A_masks_fix.mat 21.7 0 5 OP-1A_clinical_filenames.txt fix/OP-1A_clinical_fix.mat
./printBloodVelocities.sh SA SA_masks_fix.mat 25.7 0 5 SA_clinical_filenames.txt fix/SA_clinical_fix.mat
./printBloodVelocities.sh RK RK_masks_fix.mat 32.6 0 5 RK_clinical_filenames.txt fix/RK_clinical_fix.mat
./printBloodVelocities.sh SE SE_masks_fix.mat 26.0 0 5 SE_clinical_filenames.txt fix/SE_clinical_fix.mat
./printBloodVelocities.sh UB UB_masks_fix.mat 36.1 0 5 UB_clinical_filenames.txt fix/UB_clinical_fix.mat
./printBloodVelocities.sh ML ML_masks_fix.mat 95.4 0 5 ML_clinical_filenames.txt fix/ML_clinical_fix.mat
./printBloodVelocities.sh CR CR_masks_fix.mat 27.9 0 5 CR_clinical_filenames.txt fix/CR_clinical_fix.mat
./printBloodVelocities.sh LJ LJ_masks_fix.mat 17.2 1 5 LJ_clinical_filenames.txt fix/LJ_clinical_fix.mat
./printBloodVelocities.sh BM BM_masks_fix.mat 55.7 0 5 BM_KW_4_filenames.txt BM/BM_KW_4_dose.mat
