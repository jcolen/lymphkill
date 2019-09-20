#!/bin/bash
#SBATCH --time=1-00:00:00
#SBATCH --output=createMaskFiles.out
#SBATCH --error=createMaskFiles.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jc8kf@virginia.edu
#SBATCH --account=phys_nrf
#SBATCH --partition=standard
#SBATCH --mem=30000
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --tasks=1

#./createMaskTextFiles.sh BD BD_masks_short.mat BD_organs_short.mat
python createMaskTextFile.py BD

#./createMaskTextFiles.sh BF BF_masks.mat BF_organs.mat
python createMaskTextFile.py BF

#./createMaskTextFiles.sh SJ-1A SJ-1A_masks.mat SJ-1A_organs.mat
python createMaskTextFile.py SJ-1A

#./createMaskTextFiles.sh JL1 JL1_masks.mat JL1_organs.mat
python createMaskTextFile.py JL1

#./createMaskTextFiles.sh SR_1A SR_masks_short.mat SR_organs_short.mat
python createMaskTextFile.py SR_1A

#./createMaskTextFiles.sh BB BB_masks_2.mat BB_organs_new.mat
python createMaskTextFile.py BB

#./createMaskTextFiles.sh BC BC_masks.mat BC_organs.mat
python createMaskTextFile.py BC

#./createMaskTextFiles.sh HR HR_masks_short.mat HR_organs_short.mat
python createMaskTextFile.py HR

#./createMaskTextFiles.sh BM BM_masks_short.mat BM_organs_short.mat
python createMaskTextFile.py BM

#./createMaskTextFiles.sh MM1 MM1_masks.mat MM1_organs.mat
python createMaskTextFile.py MM1

#./createMaskTextFiles.sh WJ WJ_masks_short.mat WJ_organs_short.mat
python createMaskTextFile.py WJ

#./createMaskTextFiles.sh BJ BJ_masks_short.mat BJ_organs_short.mat
python createMaskTextFile.py BJ

#./createMaskTextFiles.sh WH WH_masks.mat WH_organs.mat
python createMaskTextFile.py WH

#./createMaskTextFiles.sh JL JL_masks_short.mat JL_organs_short.mat
python createMaskTextFile.py JL

#./createMaskTextFiles.sh KE-1A KE-1A_masks_short.mat KE-1A_organs_short.mat
python createMaskTextFile.py KE-1A

#./createMaskTextFiles.sh HI HI_masks_short.mat HI_organs_short.mat
python createMaskTextFile.py HI

#./createMaskTextFiles.sh ME1 ME1_masks.mat ME1_organs.mat
python createMaskTextFile.py ME1

#./createMaskTextFiles.sh ST ST_masks.mat ST_organs.mat
python createMaskTextFile.py ST

#./createMaskTextFiles.sh BR BR_masks_short.mat BR_organs_short.mat
python createMaskTextFile.py BR

#./createMaskTextFiles.sh LV-1A LV-1A_masks.mat LV-1A_organs.mat
python createMaskTextFile.py LV-1A

#./createMaskTextFiles.sh HB HB_masks.mat HB_organs.mat
python createMaskTextFile.py HB

#./createMaskTextFiles.sh RD-1A RD-1A_masks.mat RD-1A_organs.mat
python createMaskTextFile.py RD-1A

#./createMaskTextFiles.sh JF JF_masks_short.mat JF_organs_short.mat
python createMaskTextFile.py JF

#./createMaskTextFiles.sh HA HA_masks_short.mat HA_organs_short.mat
python createMaskTextFile.py HA

#./createMaskTextFiles.sh FJ FJ_masks.mat FJ_organs.mat
python createMaskTextFile.py FJ

#./createMaskTextFiles.sh HD HD_masks.mat HD_organs.mat
python createMaskTextFile.py HD

#./createMaskTextFiles.sh LD LD_masks.mat LD_organs.mat
python createMaskTextFile.py LD

#./createMaskTextFiles.sh SJ2 SJ2_masks.mat SJ2_organs.mat
python createMaskTextFile.py SJ2

#./createMaskTextFiles.sh SR1 SR1_masks.mat SR1_organs.mat
python createMaskTextFile.py SR1

#./createMaskTextFiles.sh LL LL_masks_short.mat LL_organs_short.mat
python createMaskTextFile.py LL

#./createMaskTextFiles.sh SB-1A SB-1A_masks.mat SB-1A_organs.mat
python createMaskTextFile.py SB-1A

#./createMaskTextFiles.sh HR1 HR1_masks.mat HR1_organs.mat
python createMaskTextFile.py HR1

#./createMaskTextFiles.sh BL BL_masks_short.mat BL_organs_short.mat
python createMaskTextFile.py BL

#./createMaskTextFiles.sh SH SH_masks.mat SH_organs.mat
python createMaskTextFile.py SH

#./createMaskTextFiles.sh BD1 BD1_masks_short.mat BD1_organs_short.mat
python createMaskTextFile.py BD1

#./createMaskTextFiles.sh MJ1 MJ1_masks_short.mat MJ1_organs_short.mat
python createMaskTextFile.py MJ1

#./createMaskTextFiles.sh PJ PJ_masks.mat PJ_organs.mat
python createMaskTextFile.py PJ

#./createMaskTextFiles.sh RM RM_masks.mat RM_organs.mat
python createMaskTextFile.py RM

#./createMaskTextFiles.sh PH PH_masks.mat PH_organs.mat
python createMaskTextFile.py PH

#./createMaskTextFiles.sh TC TC_masks.mat TC_organs.mat
python createMaskTextFile.py TC

#./createMaskTextFiles.sh VJ VJ_masks.mat VJ_organs.mat
python createMaskTextFile.py VJ

#./createMaskTextFiles.sh BB1 BB1_masks.mat BB1_organs.mat
python createMaskTextFile.py BB1

#./createMaskTextFiles.sh RA RA_masks.mat RA_organs.mat
python createMaskTextFile.py RA

#./createMaskTextFiles.sh AA AA_masks.mat AA_organs.mat
python createMaskTextFile.py AA

#./createMaskTextFiles.sh GV GV_masks_short.mat GV_organs_short.mat
python createMaskTextFile.py GV

#./createMaskTextFiles.sh AP AP_masks.mat AP_organs.mat
python createMaskTextFile.py AP

#./createMaskTextFiles.sh GS GS_masks_short.mat GS_organs_short.mat
python createMaskTextFile.py GS

#./createMaskTextFiles.sh CM CM_masks_short.mat CM_organs_short.mat
python createMaskTextFile.py CM

#./createMaskTextFiles.sh CL_1A CL_1A_masks.mat CL_1A_organs.mat
python createMaskTextFile.py CL_1A

#./createMaskTextFiles.sh JS-1B JS-1B_masks_short.mat JS-1B_organs_short.mat
python createMaskTextFile.py JS-1B

#./createMaskTextFiles.sh FL-1A FL-1A_masks.mat FL-1A_organs.mat
python createMaskTextFile.py FL-1A

#./createMaskTextFiles.sh CE CE_masks_short.mat CE_organs_short.mat
python createMaskTextFile.py CE

#./createMaskTextFiles.sh GF GF_masks.mat GF_organs.mat
python createMaskTextFile.py GF

#./createMaskTextFiles.sh GE GE_masks.mat GE_organs.mat
python createMaskTextFile.py GE

#./createMaskTextFiles.sh HH-1BResim HH-1BResim_masks.mat HH-1BResim_organs.mat
python createMaskTextFile.py HH-1BResim

#./createMaskTextFiles.sh EC EC_masks_short.mat EC_organs_short.mat
python createMaskTextFile.py EC

#./createMaskTextFiles.sh LR_1A LR_1A_masks_short.mat LR_1A_organs_short.mat
python createMaskTextFile.py LR_1A

#./createMaskTextFiles.sh CS CS_masks.mat CS_organs.mat
python createMaskTextFile.py CS

#./createMaskTextFiles.sh GJ GJ_masks.mat GJ_organs.mat
python createMaskTextFile.py GJ

#./createMaskTextFiles.sh WC WC_masks.mat WC_organs.mat
python createMaskTextFile.py WC

#./createMaskTextFiles.sh WJ1 WJ1_masks.mat WJ1_organs.mat
python createMaskTextFile.py WJ1

#./createMaskTextFiles.sh ME ME_masks.mat ME_organs.mat
python createMaskTextFile.py ME

#./createMaskTextFiles.sh HJ HJ_masks.mat HJ_organs_short.mat
python createMaskTextFile.py HJ

#./createMaskTextFiles.sh JM-1A JM-1A_masks_short.mat JM-1A_organs_short.mat
python createMaskTextFile.py JM-1A

#./createMaskTextFiles.sh RH RH_masks.mat RH_organs.mat
python createMaskTextFile.py RH

#./createMaskTextFiles.sh OH OH_masks.mat OH_organs.mat
python createMaskTextFile.py OH

#./createMaskTextFiles.sh MM MM_masks_short.mat MM_organs_short.mat
python createMaskTextFile.py MM

#./createMaskTextFiles.sh CA-2ALUL CA-2ALUL_masks.mat CA-2ALUL_organs.mat
python createMaskTextFile.py CA-2ALUL

#./createMaskTextFiles.sh MB MB_masks.mat MB_organs.mat
python createMaskTextFile.py MB

#./createMaskTextFiles.sh HM-1A HM-1A_masks.mat HM-1A_organs.mat
python createMaskTextFile.py HM-1A

#./createMaskTextFiles.sh SJ1-1A SJ1-1A_masks.mat SJ1-1A_organs.mat
python createMaskTextFile.py SJ1-1A

#./createMaskTextFiles.sh MK MK_masks_short.mat MK_organs_short.mat
python createMaskTextFile.py MK

#./createMaskTextFiles.sh MJ MJ_masks.mat MJ_organs.mat
python createMaskTextFile.py MJ

#./createMaskTextFiles.sh CE1 CE1_masks.mat CE1_organs.mat
python createMaskTextFile.py CE1

#./createMaskTextFiles.sh MD-1A MD-1A_masks.mat MD-1A_organs.mat
python createMaskTextFile.py MD-1A

#./createMaskTextFiles.sh MR MR_masks.mat MR_organs.mat
python createMaskTextFile.py MR

#./createMaskTextFiles.sh CJ CJ_masks.mat CJ_organs.mat
python createMaskTextFile.py CJ

#./createMaskTextFiles.sh OP-1A OP-1A_masks.mat OP-1A_organs.mat
python createMaskTextFile.py OP-1A

#./createMaskTextFiles.sh SA SA_masks.mat SA_organs.mat
python createMaskTextFile.py SA

#./createMaskTextFiles.sh RK RK_masks.mat RK_organs.mat
python createMaskTextFile.py RK

#./createMaskTextFiles.sh SE SE_masks.mat SE_organs.mat
python createMaskTextFile.py SE

#./createMaskTextFiles.sh UB UB_masks.mat UB_organs.mat
python createMaskTextFile.py UB

