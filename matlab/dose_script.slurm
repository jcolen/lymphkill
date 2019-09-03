#!/bin/bash
#SBATCH --time=1-00:00:00
#SBATCH --output=dose_script.out
#SBATCH --error=dose_script.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jc8kf@virginia.edu
#SBATCH --account=phys_nrf
#SBATCH --partition=standard

#SBATCH --mem=30000

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1

#SBATCH --tasks=1

#./dose_script.sh AR/AR_masks.mat AR/AR_clinical_filenames.txt 30.1 0 AR/AR_clinical_dose.mat
#./dose_script.sh BD/BD_masks.mat BD/BD_clinical_filenames.txt 29.3 0 BD/BD_clinical_dose.mat
#./dose_script.sh BD1/BD1_masks.mat BD1/BD1_clinical_filenames.txt 38.4 0 BD1/BD1_clinical_dose.mat
#./dose_script.sh MS/MS_masks.mat MS/MS_clinical_filenames.txt 47.1 0 MS/MS_clinical_dose.mat
#./dose_script.sh BJ/BJ_masks.mat BJ/BJ_clinical_filenames.txt 21.9 0 BJ/BJ_clinical_dose.mat
#./dose_script.sh BC/BC_masks.mat BC/BC_clinical_filenames.txt 20.0 0 BC/BC_clinical_dose.mat
#./dose_script.sh HA/HA_masks_short.mat HA/HA_clinical_filenames.txt 74.0 0 HA/HA_clinical_dose_short.mat
#./dose_script.sh BL/BL_masks.mat BL/BL_clinical_filenames.txt 62.5 0 BL/BL_clinical_dose.mat
#./dose_script.sh BM/BM_masks_short.mat BM/BM_clinical_filenames.txt 24.7 1 BM/BM_clinical_dose_short.mat
#./dose_script.sh CE/CE_masks_short.mat CE/CE_clinical_filenames.txt 38.6 0 CE/CE_clinical_dose_short.mat
#./dose_script.sh JL/JL_masks_short.mat JL/JL_clinical_filenames.txt 22.5 1 JL/JL_clinical_dose_short.mat
#./dose_script.sh AA/AA_masks.mat AA/AA_clinical_filenames.txt 46.3 0 AA/AA_clinical_dose.mat
#./dose_script.sh EC/EC_masks_short.mat EC/EC_clinical_filenames.txt 64.9 0 EC/EC_clinical_dose_short.mat
#./dose_script.sh FJ/FJ_masks.mat FJ/FJ_clinical_filenames.txt 16.9 0 FJ/FJ_clinical_dose.mat
#./dose_script.sh GV/GV_masks_short.mat GV/GV_clinical_filenames.txt 26.7 1 GV/GV_clinical_dose_short.mat
#./dose_script.sh HR1/HR1_masks.mat HR1/HR1_clinical_filenames.txt 32.0 0 HR1/HR1_clinical_dose.mat
#./dose_script.sh JF/JF_masks_short.mat JF/JF_clinical_filenames.txt 28.1 0 JF/JF_clinical_dose_short.mat
#./dose_script.sh JS-1B/JS-1B_masks_short.mat JS-1B/JS-1B_clinical_filenames.txt 22.6 0 JS-1B/JS-1B_clinical_dose_short.mat
#./dose_script.sh JL1/JL1_masks.mat JL1/JL1_clinical_filenames.txt 50.9 0 JL1/JL1_clinical_dose.mat
#./dose_script.sh KE-1A/KE-1A_masks_short.mat KE-1A/KE-1A_clinical_filenames.txt 42.6 0 KE-1A/KE-1A_clinical_dose_short.mat
#./dose_script.sh LL/LL_masks_short.mat LL/LL_clinical_filenames.txt 41.2 0 LL/LL_clinical_dose_short.mat
#./dose_script.sh LJ/LJ_masks_short.mat LJ/LJ_clinical_filenames.txt 17.2 1 LJ/LJ_clinical_dose_short.mat
#./dose_script.sh MM/MM_masks_short.mat MM/MM_clinical_filenames.txt 49.3 0 MM/MM_clinical_dose_short.mat
#./dose_script.sh JM-1A/JM-1A_masks.mat JM-1A/JM-1A_clinical_filenames.txt 38.3 1 JM-1A/JM-1A_clinical_dose.mat
#./dose_script.sh MM1/MM1_masks.mat MM1/MM1_clinical_filenames.txt 31.7 0 MM1/MM1_clinical_dose.mat
#./dose_script.sh MK/MK_masks_short.mat MK/MK_clinical_filenames.txt 18.6 0 MK/MK_clinical_dose_short.mat
#./dose_script.sh MJ1/MJ1_masks_short.mat MJ1/MJ1_clinical_filenames.txt 17.9 1 MJ1/MJ1_clinical_dose_short.mat
#./dose_script.sh GE/GE_masks.mat GE/GE_clinical_filenames.txt 40.1 0 GE/GE_clinical_dose.mat
#./dose_script.sh GF/GF_masks.mat GF/GF_clinical_filenames.txt 35.0 0 GF/GF_clinical_dose.mat
#./dose_script.sh MD-1A/MD-1A_masks.mat MD-1A/MD-1A_clinical_filenames.txt 31.9 0 MD-1A/MD-1A_clinical_dose.mat
#./dose_script.sh BB/BB_masks_2.mat BB/BB_clinical_filenames.txt 76.4 0 BB/BB_clinical_dose_2.mat
#./dose_script.sh SR_1A/SR_masks_2.mat SR_1A/SR_clinical_filenames.txt 26.7 0 SR_1A/SR_clinical_dose_2.mat
#./dose_script.sh PH/PH_masks.mat PH/PH_clinical_filenames.txt 61.1 0 PH/PH_clinical_dose.mat
#./dose_script.sh OP-1A/OP-1A_masks.mat OP-1A/OP-1A_clinical_filenames.txt 21.7 0 OP-1A/OP-1A_clinical_dose.mat
#./dose_script.sh SH/SH_masks.mat SH/SH_clinical_filenames.txt 18.5 0 SH/SH_clinical_dose.mat
#./dose_script.sh WJ1/WJ1_masks.mat WJ1/WJ1_clinical_filenames.txt 40.3 0 WJ1/WJ1_clinical_dose.mat
#./dose_script.sh WH/WH_masks.mat WH/WH_clinical_filenames.txt 25.6 0 WH/WH_clinical_dose.mat
#./dose_script.sh UB/UB_masks.mat UB/UB_clinical_filenames.txt 36.1 0 UB/UB_clinical_dose.mat
#./dose_script.sh VJ/VJ_masks.mat VJ/VJ_clinical_filenames.txt 23.2 0 VJ/VJ_clinical_dose.mat
#./dose_script.sh BB1/BB1_masks.mat BB1/BB1_clinical_filenames.txt 51.1 0 BB1/BB1_clinical_dose.mat
#./dose_script.sh ME/ME_masks.mat ME/ME_clinical_filenames.txt 38.8 0 ME/ME_clinical_dose.mat
#./dose_script.sh CA-2ALUL/CA-2ALUL_masks.mat CA-2ALUL/CA-2ALUL_clinical_filenames.txt 58.3 0 CA-2ALUL/CA-2ALUL_clinical_dose.mat
#./dose_script.sh BR/BR_masks_short.mat BR/BR_clinical_filenames.txt 24.9 0 BR/BR_clinical_dose_short.mat
#./dose_script.sh CM/CM_masks_short.mat CM/CM_clinical_filenames.txt 32.6 1 CM/CM_clinical_dose_short.mat
#./dose_script.sh GS/GS_masks_short.mat GS/GS_clinical_filenames.txt 26.9 1 GS/GS_clinical_dose_short.mat
#./dose_script.sh HR/HR_masks_short.mat HR/HR_clinical_filenames.txt 22.7 1 HR/HR_clinical_dose_short.mat
#./dose_script.sh LR_1A/LR_1A_masks_short.mat LR_1A/LR_1A_clinical_filenames.txt 28.8 0 LR_1A/LR_1A_clinical_dose_short.mat
#./dose_script.sh SR_1A/SR_masks_short.mat SR_1A/SR_clinical_filenames.txt 26.7 0 SR_1A/SR_clinical_dose_short.mat
#./dose_script.sh WJ/WJ_masks_short.mat WJ/WJ_clinical_filenames.txt 47.9 0 WJ/WJ_clinical_dose_short.mat
#./dose_script.sh HI/HI_masks_short.mat HI/HI_clinical_filenames.txt 136.9 0 HI/HI_clinical_dose_short.mat
#./dose_script.sh BD/BD_masks_short.mat BD/BD_clinical_filenames.txt 29.3 0 BD/BD_clinical_dose_short.mat
#./dose_script.sh BD1/BD1_masks_short.mat BD1/BD1_clinical_filenames.txt 38.4 0 BD1/BD1_clinical_dose_short.mat
#./dose_script.sh BJ/BJ_masks_short.mat BJ/BJ_clinical_filenames.txt 21.9 0 BJ/BJ_clinical_dose_short.mat
#./dose_script.sh BL/BL_masks_short.mat BL/BL_clinical_filenames.txt 62.5 0 BL/BL_clinical_dose_short.mat
#./dose_script.sh CL_1A/CL_1A_masks.mat CL_1A/CL_1A_clinical_filenames.txt 25.2 0 CL_1A/CL_1A_clinical_dose.mat
#./dose_script.sh JM-1A/JM-1A_masks_short.mat JM-1A/JM-1A_clinical_filenames.txt 38.3 1 JM-1A/JM-1A_clinical_dose_short.mat
#./dose_script.sh OH/OH_masks.mat OH/OH_clinical_filenames.txt 49.3 0 OH/OH_clinical_dose.mat
#./dose_script.sh RA/RA_masks.mat RA/RA_clinical_filenames.txt 21.5 0 RA/RA_clinical_dose.mat
#./dose_script.sh RK/RK_masks.mat RK/RK_clinical_filenames.txt 32.6 0 RK/RK_clinical_dose.mat
#./dose_script.sh SA/SA_masks.mat SA/SA_clinical_filenames.txt 25.7 0 SA/SA_clinical_dose.mat
#./dose_script.sh SJ1-1A/SJ1-1A_masks.mat SJ1-1A/SJ1-1A_clinical_filenames.txt 40.5 1 SJ1-1A/SJ1-1A_clinical_dose.mat
#./dose_script.sh SJ2/SJ2_masks.mat SJ2/SJ2_clinical_filenames.txt 55.2 0 SJ2/SJ2_clinical_dose.mat
#./dose_script.sh TC/TC_masks.mat TC/TC_clinical_filenames.txt 26.2 0 TC/TC_clinical_dose.mat
#./dose_script.sh WC/WC_masks.mat WC/WC_clinical_filenames.txt 20.0 0 WC/WC_clinical_dose.mat
#./dose_script.sh SR1/SR1_masks.mat SR1/SR1_clinical_filenames.txt 53.6 0 SR1/SR1_clinical_dose.mat
#./dose_script.sh MJ/MJ_masks.mat MJ/MJ_clinical_filenames.txt 25.2 0 MJ/MJ_clinical_dose.mat
#./dose_script.sh HJ/HJ_masks.mat HJ/HJ_clinical_filenames.txt 53.2 0 HJ/HJ_clinical_dose.mat
#./dose_script.sh RD-1A/RD-1A_masks.mat RD-1A/RD-1A_clinical_filenames.txt 37.7 0 RD-1A/RD-1A_clinical_dose.mat
#./dose_script.sh RM/RM_masks.mat RM/RM_clinical_filenames.txt 34.3 0 RM/RM_clinical_dose.mat
#./dose_script.sh PJ/PJ_masks.mat PJ/PJ_clinical_filenames.txt 14.3 0 PJ/PJ_clinical_dose.mat
#./dose_script.sh AP/AP_masks.mat AP/AP_clinical_filenames.tx 28.4 0 AP/AP_clinical_dose.mat
#./dose_script.sh BF/BF_masks.mat BF/BF_clinical_filenames.tx 29.0 0 BF/BF_clinical_dose.mat
#./dose_script.sh CE1/CE1_masks.mat CE1/CE1_clinical_filenames.tx 28.2 0 CE1/CE1_clinical_dose.mat
#./dose_script.sh CJ/CJ_masks.mat CJ/CJ_clinical_filenames.tx 31.4 1 CJ/CJ_clinical_dose.mat
#./dose_script.sh CS/CS_masks.mat CS/CS_clinical_filenames.tx 66.9 0 CS/CS_clinical_dose.mat
#./dose_script.sh FL-1A/FL-1A_masks.mat FL-1A/FL-1A_clinical_filenames.tx 52.2 0 FL-1A/FL-1A_clinical_dose.mat
#./dose_script.sh GJ/GJ_masks.mat GJ/GJ_clinical_filenames.tx 35.4 0 GJ/GJ_clinical_dose.mat
#./dose_script.sh HB/HB_masks.mat HB/HB_clinical_filenames.tx 35.4 0 HB/HB_clinical_dose.mat
#./dose_script.sh HD/HD_masks.mat HD/HD_clinical_filenames.tx 24.3 0 HD/HD_clinical_dose.mat
#./dose_script.sh HH-1BResim/HH-1BResim_masks.mat HH-1BResim/HH-1BResim_clinical_filenames.tx 43.4 1 HH-1BResim/HH-1BResim_clinical_dose.mat
#./dose_script.sh HM-1A/HM-1A_masks.mat HM-1A/HM-1A_clinical_filenames.tx 29.3 1 HM-1A/HM-1A_clinical_dose.mat
#./dose_script.sh LD/LD_masks.mat LD/LD_clinical_filenames.tx 20.7 1 LD/LD_clinical_dose.mat
#./dose_script.sh LV-1A/LV-1A_masks.mat LV-1A/LV-1A_clinical_filenames.tx 34.8 0 LV-1A/LV-1A_clinical_dose.mat
#./dose_script.sh MB/MB_masks.mat MB/MB_clinical_filenames.tx 40.4 0 MB/MB_clinical_dose.mat
#./dose_script.sh ME1/ME1_masks.mat ME1/ME1_clinical_filenames.tx 29.2 1 ME1/ME1_clinical_dose.mat
#./dose_script.sh MR/MR_masks.mat MR/MR_clinical_filenames.tx 42.5 0 MR/MR_clinical_dose.mat
#./dose_script.sh RH/RH_masks.mat RH/RH_clinical_filenames.tx 56.1 0 RH/RH_clinical_dose.mat
#./dose_script.sh SB-1A/SB-1A_masks.mat SB-1A/SB-1A_clinical_filenames.tx 27.7 0 SB-1A/SB-1A_clinical_dose.mat
#./dose_script.sh SE/SE_masks.mat SE/SE_clinical_filenames.tx 26.0  0 SE/SE_clinical_dose.mat
#./dose_script.sh SJ-1A/SJ-1A_masks.mat SJ-1A/SJ-1A_clinical_filenames.tx 33.3 1 SJ-1A/SJ-1A_clinical_dose.mat
#./dose_script.sh ST/ST_masks.mat ST/ST_clinical_filenames.tx 32.6 1 ST/ST_clinical_dose.mat
#./dose_script.sh RC/RC_masks.mat RC/RC_clinical_filenames.txt 27.5 0 RC/RC_clinical_dose.mat
#./dose_script.sh RC/RC_masks.mat RC/RC_replan_filenames.txt 93.4 0 RC/RC_replan_1400_dose.mat
#./dose_script.sh RC/RC_masks.mat RC/RC_replan_filenames.txt 54.5 0 RC/RC_replan_2400_dose.mat
#./dose_script.sh RC/RC_masks.mat RC/RC_replan_filenames.txt 93.4 1 RC/RC_replan_1400_gated_dose.mat
#./dose_script.sh RC/RC_masks.mat RC/RC_replan_filenames.txt 54.5 1 RC/RC_replan_2400_gated_dose.mat
