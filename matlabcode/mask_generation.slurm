#!/bin/bash
#SBATCH --time=1-00:00:00
#SBATCH --output=mask_generation.out
#SBATCH --error=mask_generation.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jc8kf@virginia.edu
#SBATCH --account=phys_nrf
#SBATCH --partition=standard
#SBATCH --mem=30000
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --tasks=1

./mask_generation.sh BD BD_organs_short.mat BD_dicom_info.mat BD_masks_fix.mat
./mask_generation.sh BF BF_organs.mat BF_dicom_info.mat BF_masks_fix.mat
./mask_generation.sh SJ-1A SJ-1A_organs.mat SJ-1A_dicom_info.mat SJ-1A_masks_fix.mat
./mask_generation.sh JL1 JL1_organs.mat JL1_dicom_info.mat JL1_masks_fix.mat
./mask_generation.sh SR_1A SR_organs_short.mat SR_1A_dicom_info.mat SR_1A_masks_fix.mat
./mask_generation.sh BB BB_organs_new.mat BB_dicom_info.mat BB_masks_fix.mat
./mask_generation.sh BC BC_organs.mat BC_dicom_info.mat BC_masks_fix.mat
./mask_generation.sh HR HR_organs_short.mat HR_dicom_info.mat HR_masks_fix.mat
./mask_generation.sh BM BM_organs_short.mat BM_dicom_info.mat BM_masks_fix.mat
./mask_generation.sh MM1 MM1_organs.mat MM1_dicom_info.mat MM1_masks_fix.mat
./mask_generation.sh WJ WJ_organs_short.mat WJ_dicom_info.mat WJ_masks_fix.mat
./mask_generation.sh BJ BJ_organs_short.mat BJ_dicom_info.mat BJ_masks_fix.mat
./mask_generation.sh WH WH_organs.mat WH_dicom_info.mat WH_masks_fix.mat
./mask_generation.sh JL JL_organs_short.mat JL_dicom_info.mat JL_masks_fix.mat
./mask_generation.sh KE-1A KE-1A_organs_short.mat KE-1A_dicom_info.mat KE-1A_masks_fix.mat
./mask_generation.sh HI HI_organs_short.mat HI_dicom_info.mat HI_masks_fix.mat
./mask_generation.sh ME1 ME1_organs.mat ME1_dicom_info.mat ME1_masks_fix.mat
./mask_generation.sh ST ST_organs.mat ST_dicom_info.mat ST_masks_fix.mat
./mask_generation.sh BR BR_organs_short.mat BR_dicom_info.mat BR_masks_fix.mat
./mask_generation.sh LV-1A LV-1A_organs.mat LV-1A_dicom_info.mat LV-1A_masks_fix.mat
./mask_generation.sh HB HB_organs.mat HB_dicom_info.mat HB_masks_fix.mat
./mask_generation.sh RD-1A RD-1A_organs.mat RD-1A_dicom_info.mat RD-1A_masks_fix.mat
./mask_generation.sh JF JF_organs_short.mat JF_dicom_info.mat JF_masks_fix.mat
./mask_generation.sh HA HA_organs_short.mat HA_dicom_info.mat HA_masks_fix.mat
./mask_generation.sh FJ FJ_organs.mat FJ_dicom_info.mat FJ_masks_fix.mat
./mask_generation.sh HD HD_organs.mat HD_dicom_info.mat HD_masks_fix.mat
./mask_generation.sh LD LD_organs.mat LD_dicom_info.mat LD_masks_fix.mat
./mask_generation.sh SJ2 SJ2_organs.mat SJ2_dicom_info.mat SJ2_masks_fix.mat
./mask_generation.sh SR1 SR1_organs.mat SR1_dicom_info.mat SR1_masks_fix.mat
./mask_generation.sh LL LL_organs_short.mat LL_dicom_info.mat LL_masks_fix.mat
./mask_generation.sh SB-1A SB-1A_organs.mat SB-1A_dicom_info.mat SB-1A_masks_fix.mat
./mask_generation.sh HR1 HR1_organs.mat HR1_dicom_info.mat HR1_masks_fix.mat
./mask_generation.sh BL BL_organs_short.mat BL_dicom_info.mat BL_masks_fix.mat
./mask_generation.sh SH SH_organs.mat SH_dicom_info.mat SH_masks_fix.mat
./mask_generation.sh BD1 BD1_organs_short.mat BD1_dicom_info.mat BD1_masks_fix.mat
./mask_generation.sh MJ1 MJ1_organs_short.mat MJ1_dicom_info.mat MJ1_masks_fix.mat
./mask_generation.sh PJ PJ_organs.mat PJ_dicom_info.mat PJ_masks_fix.mat
./mask_generation.sh RM RM_organs.mat RM_dicom_info.mat RM_masks_fix.mat
./mask_generation.sh PH PH_organs.mat PH_dicom_info.mat PH_masks_fix.mat
./mask_generation.sh TC TC_organs.mat TC_dicom_info.mat TC_masks_fix.mat
./mask_generation.sh VJ VJ_organs.mat VJ_dicom_info.mat VJ_masks_fix.mat
./mask_generation.sh BB1 BB1_organs.mat BB1_dicom_info.mat BB1_masks_fix.mat
./mask_generation.sh RA RA_organs.mat RA_dicom_info.mat RA_masks_fix.mat
./mask_generation.sh AA AA_organs.mat AA_dicom_info.mat AA_masks_fix.mat
./mask_generation.sh GV GV_organs_short.mat GV_dicom_info.mat GV_masks_fix.mat
./mask_generation.sh AP AP_organs.mat AP_dicom_info.mat AP_masks_fix.mat
./mask_generation.sh GS GS_organs_short.mat GS_dicom_info.mat GS_masks_fix.mat
./mask_generation.sh CM CM_organs_short.mat CM_dicom_info.mat CM_masks_fix.mat
./mask_generation.sh CL_1A CL_1A_organs.mat CL_1A_dicom_info.mat CL_1A_masks_fix.mat
./mask_generation.sh JS-1B JS-1B_organs_short.mat JS-1B_dicom_info.mat JS-1B_masks_fix.mat
./mask_generation.sh FL-1A FL-1A_organs.mat FL-1A_dicom_info.mat FL-1A_masks_fix.mat
./mask_generation.sh CE CE_organs_short.mat CE_dicom_info.mat CE_masks_fix.mat
./mask_generation.sh GF GF_organs.mat GF_dicom_info.mat GF_masks_fix.mat
./mask_generation.sh GE GE_organs.mat GE_dicom_info.mat GE_masks_fix.mat
./mask_generation.sh HH-1BResim HH-1BResim_organs.mat HH-1BResim_dicom_info.mat HH-1BResim_masks_fix.mat
./mask_generation.sh EC EC_organs_short.mat EC_dicom_info.mat EC_masks_fix.mat
./mask_generation.sh LR_1A LR_1A_organs_short.mat LR_1A_dicom_info.mat LR_1A_masks_fix.mat
./mask_generation.sh CS CS_organs.mat CS_dicom_info.mat CS_masks_fix.mat
./mask_generation.sh GJ GJ_organs.mat GJ_dicom_info.mat GJ_masks_fix.mat
./mask_generation.sh WC WC_organs.mat WC_dicom_info.mat WC_masks_fix.mat
./mask_generation.sh WJ1 WJ1_organs.mat WJ1_dicom_info.mat WJ1_masks_fix.mat
./mask_generation.sh ME ME_organs.mat ME_dicom_info.mat ME_masks_fix.mat
./mask_generation.sh HJ HJ_organs_short.mat HJ_dicom_info.mat HJ_masks_fix.mat
./mask_generation.sh JM-1A JM-1A_organs_short.mat JM-1A_dicom_info.mat JM-1A_masks_fix.mat
./mask_generation.sh RH RH_organs.mat RH_dicom_info.mat RH_masks_fix.mat
./mask_generation.sh OH OH_organs.mat OH_dicom_info.mat OH_masks_fix.mat
./mask_generation.sh MM MM_organs_short.mat MM_dicom_info.mat MM_masks_fix.mat
./mask_generation.sh CA-2ALUL CA-2ALUL_organs.mat CA-2ALUL_dicom_info.mat CA-2ALUL_masks_fix.mat
./mask_generation.sh MB MB_organs.mat MB_dicom_info.mat MB_masks_fix.mat
./mask_generation.sh HM-1A HM-1A_organs.mat HM-1A_dicom_info.mat HM-1A_masks_fix.mat
./mask_generation.sh SJ1-1A SJ1-1A_organs.mat SJ1-1A_dicom_info.mat SJ1-1A_masks_fix.mat
./mask_generation.sh MK MK_organs_short.mat MK_dicom_info.mat MK_masks_fix.mat
./mask_generation.sh MJ MJ_organs.mat MJ_dicom_info.mat MJ_masks_fix.mat
./mask_generation.sh CE1 CE1_organs.mat CE1_dicom_info.mat CE1_masks_fix.mat
./mask_generation.sh MD-1A MD-1A_organs.mat MD-1A_dicom_info.mat MD-1A_masks_fix.mat
./mask_generation.sh MR MR_organs.mat MR_dicom_info.mat MR_masks_fix.mat
./mask_generation.sh CJ CJ_organs.mat CJ_dicom_info.mat CJ_masks_fix.mat
./mask_generation.sh OP-1A OP-1A_organs.mat OP-1A_dicom_info.mat OP-1A_masks_fix.mat
./mask_generation.sh SA SA_organs.mat SA_dicom_info.mat SA_masks_fix.mat
./mask_generation.sh RK RK_organs.mat RK_dicom_info.mat RK_masks_fix.mat
./mask_generation.sh SE SE_organs.mat SE_dicom_info.mat SE_masks_fix.mat
./mask_generation.sh UB UB_organs.mat UB_dicom_info.mat UB_masks_fix.mat
./mask_generation.sh ML ML_organs_new.mat ML_dicom_info.mat ML_masks_fix.mat
./mask_generation.sh CR CR_organs.mat CR_dicom_info.mat CR_masks_fix.mat
./mask_generation.sh MS MS_organs.mat MS_dicom_info.mat MS_masks_fix.mat
./mask_generation.sh LJ LJ_organs.mat LJ_dicom_info.mat LJ_masks_fix.mat
./mask_generation.sh AR AR_organs.mat AR_dicom_info.mat AR_masks_fix.mat
