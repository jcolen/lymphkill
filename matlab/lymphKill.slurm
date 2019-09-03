#!/bin/bash
#SBATCH --time=5-00:00:00
#SBATCH --output=lymphKill.out
#SBATCH --error=lymphKill.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jc8kf@virginia.edu
#SBATCH --account=phys_nrf
#SBATCH --partition=standard
#SBATCH --mem=30000
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --tasks=1

./lymphKill.sh BD BD_masks_fix.mat 29.3 0 5 BD_clinical_filenames.txt fix/BD_clinical_fix.mat
./lymphKill.sh BF BF_masks_fix.mat 29.0 0 5 BF_clinical_filenames.tx fix/BF_clinical_fix.mat
./lymphKill.sh SJ-1A SJ-1A_masks_fix.mat 33.3 1 5 SJ-1A_clinical_filenames.tx fix/SJ-1A_clinical_fix.mat
./lymphKill.sh JL1 JL1_masks_fix.mat 50.9 0 5 JL1_clinical_filenames.txt fix/JL1_clinical_fix.mat
./lymphKill.sh SR_1A SR_1A_masks_fix.mat 26.7 0 5 SR_clinical_filenames.txt fix/SR_1A_clinical_fix.mat
./lymphKill.sh BB BB_masks_fix.mat 76.4 0 5 BB_clinical_filenames.txt fix/BB_clinical_fix.mat
./lymphKill.sh BC BC_masks_fix.mat 20.0 0 5 BC_clinical_filenames.txt fix/BC_clinical_fix.mat
./lymphKill.sh HR HR_masks_fix.mat 22.7 1 5 HR_clinical_filenames.txt fix/HR_clinical_fix.mat
./lymphKill.sh BM BM_masks_fix.mat 24.7 1 5 BM_clinical_filenames.txt fix/BM_clinical_fix.mat
./lymphKill.sh MM1 MM1_masks_fix.mat 31.7 0 5 MM1_clinical_filenames.txt fix/MM1_clinical_fix.mat
./lymphKill.sh WJ WJ_masks_fix.mat 47.9 0 5 WJ_clinical_filenames.txt fix/WJ_clinical_fix.mat
./lymphKill.sh BJ BJ_masks_fix.mat 21.9 0 5 BJ_clinical_filenames.txt fix/BJ_clinical_fix.mat
./lymphKill.sh WH WH_masks_fix.mat 25.6 0 5 WH_clinical_filenames.txt fix/WH_clinical_fix.mat
./lymphKill.sh JL JL_masks_fix.mat 22.5 1 5 JL_clinical_filenames.txt fix/JL_clinical_fix.mat
./lymphKill.sh KE-1A KE-1A_masks_fix.mat 42.6 0 5 KE-1A_clinical_filenames.txt fix/KE-1A_clinical_fix.mat
./lymphKill.sh HI HI_masks_fix.mat 136.9 0 5 HI_clinical_filenames.txt fix/HI_clinical_fix.mat
./lymphKill.sh ME1 ME1_masks_fix.mat 29.2 1 5 ME1_clinical_filenames.tx fix/ME1_clinical_fix.mat
./lymphKill.sh HH-1BResim HH-1BResim_masks_fix.mat 43.4 1 5 HH-1BResim_clinical_filenames.tx fix/HH-1BResim_clinical_fix.mat
./lymphKill.sh BR BR_masks_fix.mat 24.9 0 5 BR_clinical_filenames.txt fix/BR_clinical_fix.mat
./lymphKill.sh LV-1A LV-1A_masks_fix.mat 34.8 0 5 LV-1A_clinical_filenames.tx fix/LV-1A_clinical_fix.mat
./lymphKill.sh HB HB_masks_fix.mat 35.4 0 5 HB_clinical_filenames.tx fix/HB_clinical_fix.mat
./lymphKill.sh WC WC_masks_fix.mat 20.0 0 5 WC_clinical_filenames.txt fix/WC_clinical_fix.mat
./lymphKill.sh FJ FJ_masks_fix.mat 16.9 0 5 FJ_clinical_filenames.txt fix/FJ_clinical_fix.mat
./lymphKill.sh HA HA_masks_fix.mat 74.0 0 5 HA_clinical_filenames.txt fix/HA_clinical_fix.mat
./lymphKill.sh JF JF_masks_fix.mat 28.1 0 5 JF_clinical_filenames.txt fix/JF_clinical_fix.mat
./lymphKill.sh HD HD_masks_fix.mat 24.3 0 5 HD_clinical_filenames.tx fix/HD_clinical_fix.mat
./lymphKill.sh LD LD_masks_fix.mat 20.7 1 5 LD_clinical_filenames.tx fix/LD_clinical_fix.mat
./lymphKill.sh SJ2 SJ2_masks_fix.mat 55.2 0 5 SJ2_clinical_filenames.txt fix/SJ2_clinical_fix.mat
./lymphKill.sh SR1 SR1_masks_fix.mat 53.6 0 5 SR1_clinical_filenames.txt fix/SR1_clinical_fix.mat
./lymphKill.sh LL LL_masks_fix.mat 41.2 0 5 LL_clinical_filenames.txt fix/LL_clinical_fix.mat
./lymphKill.sh SB-1A SB-1A_masks_fix.mat 27.7 0 5 SB-1A_clinical_filenames.tx fix/SB-1A_clinical_fix.mat
./lymphKill.sh HR1 HR1_masks_fix.mat 32.0 0 5 HR1_clinical_filenames.txt fix/HR1_clinical_fix.mat
./lymphKill.sh BL BL_masks_fix.mat 62.5 0 5 BL_clinical_filenames.txt fix/BL_clinical_fix.mat
./lymphKill.sh RD-1A RD-1A_masks_fix.mat 37.7 0 5 RD-1A_clinical_filenames.txt fix/RD-1A_clinical_fix.mat
./lymphKill.sh SH SH_masks_fix.mat 18.5 0 5 SH_clinical_filenames.txt fix/SH_clinical_fix.mat
./lymphKill.sh BD1 BD1_masks_fix.mat 38.4 0 5 BD1_clinical_filenames.txt fix/BD1_clinical_fix.mat
./lymphKill.sh MJ1 MJ1_masks_fix.mat 17.9 1 5 MJ1_clinical_filenames.txt fix/MJ1_clinical_fix.mat
./lymphKill.sh PJ PJ_masks_fix.mat 14.3 0 5 PJ_clinical_filenames.txt fix/PJ_clinical_fix.mat
./lymphKill.sh RM RM_masks_fix.mat 34.3 0 5 RM_clinical_filenames.txt fix/RM_clinical_fix.mat
./lymphKill.sh PH PH_masks_fix.mat 61.1 0 5 PH_clinical_filenames.txt fix/PH_clinical_fix.mat
./lymphKill.sh RH RH_masks_fix.mat 56.1 0 5 RH_clinical_filenames.tx fix/RH_clinical_fix.mat
./lymphKill.sh VJ VJ_masks_fix.mat 23.2 0 5 VJ_clinical_filenames.txt fix/VJ_clinical_fix.mat
./lymphKill.sh BB1 BB1_masks_fix.mat 51.1 0 5 BB1_clinical_filenames.txt fix/BB1_clinical_fix.mat
./lymphKill.sh RA RA_masks_fix.mat 21.5 0 5 RA_clinical_filenames.txt fix/RA_clinical_fix.mat
./lymphKill.sh AA AA_masks_fix.mat 46.3 0 5 AA_clinical_filenames.txt fix/AA_clinical_fix.mat
./lymphKill.sh GV GV_masks_fix.mat 26.7 1 5 GV_clinical_filenames.txt fix/GV_clinical_fix.mat
./lymphKill.sh AP AP_masks_fix.mat 28.4 0 5 AP_clinical_filenames.tx fix/AP_clinical_fix.mat
./lymphKill.sh GS GS_masks_fix.mat 26.9 1 5 GS_clinical_filenames.txt fix/GS_clinical_fix.mat
./lymphKill.sh CM CM_masks_fix.mat 32.6 1 5 CM_clinical_filenames.txt fix/CM_clinical_fix.mat
./lymphKill.sh CL_1A CL_1A_masks_fix.mat 25.2 0 5 CL_1A_clinical_filenames.txt fix/CL_1A_clinical_fix.mat
./lymphKill.sh JS-1B JS-1B_masks_fix.mat 22.6 0 5 JS-1B_clinical_filenames.txt fix/JS-1B_clinical_fix.mat
./lymphKill.sh FL-1A FL-1A_masks_fix.mat 52.2 0 5 FL-1A_clinical_filenames.tx fix/FL-1A_clinical_fix.mat
./lymphKill.sh CE CE_masks_fix.mat 38.6 0 5 CE_clinical_filenames.txt fix/CE_clinical_fix.mat
./lymphKill.sh MJ MJ_masks_fix.mat 25.2 0 5 MJ_clinical_filenames.txt fix/MJ_clinical_fix.mat
./lymphKill.sh GF GF_masks_fix.mat 35.0 0 5 GF_clinical_filenames.txt fix/GF_clinical_fix.mat
./lymphKill.sh GE GE_masks_fix.mat 40.1 0 5 GE_clinical_filenames.txt fix/GE_clinical_fix.mat
./lymphKill.sh AR AR_masks_fix.mat 30.1 0 5 AR_clinical_filenames.txt fix/AR_clinical_fix.mat
./lymphKill.sh EC EC_masks_fix.mat 64.9 0 5 EC_clinical_filenames.txt fix/EC_clinical_fix.mat
./lymphKill.sh LR_1A LR_1A_masks_fix.mat 28.8 0 5 LR_1A_clinical_filenames.txt fix/LR_1A_clinical_fix.mat
./lymphKill.sh CS CS_masks_fix.mat 66.9 0 5 CS_clinical_filenames.tx fix/CS_clinical_fix.mat
./lymphKill.sh GJ GJ_masks_fix.mat 35.4 0 5 GJ_clinical_filenames.tx fix/GJ_clinical_fix.mat
./lymphKill.sh WJ1 WJ1_masks_fix.mat 40.3 0 5 WJ1_clinical_filenames.txt fix/WJ1_clinical_fix.mat
./lymphKill.sh ME ME_masks_fix.mat 38.8 0 5 ME_clinical_filenames.txt fix/ME_clinical_fix.mat
./lymphKill.sh HJ HJ_masks_fix.mat 53.2 0 5 HJ_clinical_filenames.txt fix/HJ_clinical_fix.mat
./lymphKill.sh JM-1A JM-1A_masks_fix.mat 38.3 1 3 JM-1A_clinical_filenames.txt fix/JM-1A_clinical_fix.mat
./lymphKill.sh ST ST_masks_fix.mat 32.6 1 5 ST_clinical_filenames.tx fix/ST_clinical_fix.mat
./lymphKill.sh OH OH_masks_fix.mat 49.3 0 5 OH_clinical_filenames.txt fix/OH_clinical_fix.mat
./lymphKill.sh MM MM_masks_fix.mat 49.3 0 5 MM_clinical_filenames.txt fix/MM_clinical_fix.mat
./lymphKill.sh CA-2ALUL CA-2ALUL_masks_fix.mat 58.3 0 5 CA-2ALUL_clinical_filenames.txt fix/CA-2ALUL_clinical_fix.mat
./lymphKill.sh MB MB_masks_fix.mat 40.4 0 3 MB_clinical_filenames.tx fix/MB_clinical_fix.mat
./lymphKill.sh HM-1A HM-1A_masks_fix.mat 29.3 1 5 HM-1A_clinical_filenames.tx fix/HM-1A_clinical_fix.mat
./lymphKill.sh SJ1-1A SJ1-1A_masks_fix.mat 40.5 1 5 SJ1-1A_clinical_filenames.txt fix/SJ1-1A_clinical_fix.mat
./lymphKill.sh MK MK_masks_fix.mat 18.6 0 5 MK_clinical_filenames.txt fix/MK_clinical_fix.mat
./lymphKill.sh TC TC_masks_fix.mat 26.2 0 5 TC_clinical_filenames.txt fix/TC_clinical_fix.mat
./lymphKill.sh CE1 CE1_masks_fix.mat 28.2 0 5 CE1_clinical_filenames.txt fix/CE1_clinical_fix.mat
./lymphKill.sh MD-1A MD-1A_masks_fix.mat 31.9 0 3 MD-1A_clinical_filenames.txt fix/MD-1A_clinical_fix.mat
./lymphKill.sh MS MS_masks_fix.mat 47.1 0 5 MS_clinical_filenames.txt fix/MS_clinical_fix.mat
./lymphKill.sh MR MR_masks_fix.mat 42.5 0 5 MR_clinical_filenames.tx fix/MR_clinical_fix.mat
./lymphKill.sh CJ CJ_masks_fix.mat 31.4 1 5 CJ_clinical_filenames.tx fix/CJ_clinical_fix.mat
./lymphKill.sh OP-1A OP-1A_masks_fix.mat 21.7 0 5 OP-1A_clinical_filenames.txt fix/OP-1A_clinical_fix.mat
./lymphKill.sh SA SA_masks_fix.mat 25.7 0 5 SA_clinical_filenames.txt fix/SA_clinical_fix.mat
./lymphKill.sh RK RK_masks_fix.mat 32.6 0 5 RK_clinical_filenames.txt fix/RK_clinical_fix.mat
./lymphKill.sh SE SE_masks_fix.mat 26.0 0 5 SE_clinical_filenames.txt fix/SE_clinical_fix.mat
./lymphKill.sh UB UB_masks_fix.mat 36.1 0 5 UB_clinical_filenames.txt fix/UB_clinical_fix.mat
./lymphKill.sh ML ML_masks_fix.mat 95.4 0 5 ML_clinical_filenames.txt fix/ML_clinical_fix.mat
./lymphKill.sh CR CR_masks_fix.mat 27.9 0 5 CR_clinical_filenames.txt fix/CR_clinical_fix.mat
./lymphKill.sh LJ LJ_masks_fix.mat 17.2 1 5 LJ_clinical_filenames.txt fix/LJ_clinical_fix.mat
#./lymphKill.sh BM BM_masks_fix.mat 55.7 0 5 BM_KW_4_filenames.txt BM/BM_KW_4_dose.mat
#./lymphKill.sh BM BM_masks_fix.mat 60.0 0 4 BM_proton_RP_RD_filenames.txt BM/BM_proton_RP_RD_dose.mat --fix
#./lymphKill.sh CR CR_masks_fix.mat 85.4 0 5 CR_KW2_filenames.txt CR/CR_KW2_dose.mat --fix
#./lymphKill.sh LL LL_masks_fix.mat 45.8 0 5 LL_4Arc_EA_filenames.txt LL/LL_4Arc_EA_dose.mat --fix
#./lymphKill.sh ML ML_masks_fix.mat 99.8 0 5 ML_Replan_2arcs_filenames.txt ML/ML_Replan_2arcs_dose.mat
#./lymphKill.sh ML ML_masks_fix.mat 60.0 0 4 ML_proton_RP_RD_filenames.txt ML/ML_proton_RP_RD_dose.mat --fix
#./lymphKill.sh RC RC_masks_fix.mat 45.4 0 5 RC_KW_FullArc24_filenames.txt RC/RC_KW_FullArc24_dose.mat --fix

