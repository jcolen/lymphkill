#!/bin/bash
#SBATCH --time=5-00:00:00
#SBATCH --output=run_pipeline.out
#SBATCH --error=run_pipeline.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kn6vv@virginia.edu
#SBATCH --account=phys_nrf
#SBATCH --partition=standard
#SBATCH --mem=30000
#SBATCH --nodes=1
#SBATCH --exclusive

module unload python
module load gcc/7.1.0
module load intel/18.0
module load mvapich2/2.3.1
module load openmpi/3.1.4
module load python/3.6.6

source /nv/vol141/phys_nrf/kevin/proj/env/bin/activate


#python lymphkill/run_pipeline.py ProductionPatients/BE/
#python lymphill/run_pipeline.py NonLungSBRT/ 
#python lymphkill/run_pipeline.py data/PH/ 
#python lymphkill/run_pipeline.py data/RM/ 
#python lymphkill/run_pipeline.py data/ST/ 
#python lymphkill/run_pipeline.py data/TC/ 
#python lymphkill/run_pipeline.py data/WC/ 
#python lymphkill/run_pipeline.py data/HI/ 


#echo patient LL
#python lymphkill/structure_loading.py Cardiac/LL/
#python lymphkill/mask_generation.py Cardiac/LL/
#python get_voxel_doses.py Cardiac/LL/

#echo patient BE 
#python lymphkill/structure_loading.py Cardiac/BE/
#python lymphkill/mask_generation.py Cardiac/BE/
#python get_voxel_doses.py Cardiac/BE/


#echo patient CW
#python lymphkill/structure_loading.py Cardiac/CW/
#python lymphkill/mask_generation.py Cardiac/CW/
#python get_voxel_doses.py Cardiac/CW/


#echo patient DJ1
#python lymphkill/structure_loading.py Cardiac/DJ1/
#python lymphkill/mask_generation.py Cardiac/DJ1/
#python get_voxel_doses.py Cardiac/DJ1/


#echo patient HM1
#python lymphkill/structure_loading.py Cardiac/HM1/
#python lymphkill/mask_generation.py Cardiac/HM1/
#python get_voxel_doses.py Cardiac/HM1/


#echo patient LJ
#python lymphkill/structure_loading.py Cardiac/LJ/
#python lymphkill/mask_generation.py Cardiac/LJ/
#python get_voxel_doses.py Cardiac/LJ/


#echo patient TC
#python lymphkill/structure_loading.py Cardiac/TC/
#python lymphkill/mask_generation.py Cardiac/TC/
#python get_voxel_doses.py Cardiac/TC/





#Prodcution PAtients

#echo patient CW
#python lymphkill/structure_loading.py ProductionPatients/CW/
#python lymphkill/mask_generation.py ProductionPatients/CW/
#python get_voxel_doses.py ProductionPatients/CW/

#echo patient DJ1
#python lymphkill/structure_loading.py ProductionPatients/DJ1/
#python lymphkill/mask_generation.py ProductionPatients/DJ1/
#python get_voxel_doses.py ProductionPatients/DJ1/

#echo patient HM1
#python lymphkill/structure_loading.py ProductionPatients/HM1/
#python lymphkill/mask_generation.py ProductionPatients/HM1/
#python get_voxel_doses.py ProductionPatients/HM1/

#echo patient LJ
#python lymphkill/structure_loading.py ProductionPatients/LJ/
#python lymphkill/mask_generation.py ProductionPatients/LJ/
#python get_voxel_doses.py ProductionPatients/LJ/

#echo patient MM1
#python lymphkill/structure_loading.py ProductionPatients/MM1/
#python lymphkill/mask_generation.py ProductionPatients/MM1/
#python get_voxel_doses.py ProductionPatients/MM1/

#echo patient PJ2
#python lymphkill/structure_loading.py ProductionPatients/PJ2/
#python lymphkill/mask_generation.py ProductionPatients/PJ2/
#python get_voxel_doses.py ProductionPatients/PJ2/

#echo patient WJ
#python lymphkill/structure_loading.py ProductionPatients/WJ/
#python lymphkill/mask_generation.py ProductionPatients/WJ/
#python get_voxel_doses.py ProductionPatients/WJ/

#echo patient WL
#python lymphkill/structure_loading.py ProductionPatients/WL/
#python lymphkill/mask_generation.py ProductionPatients/WL/
#python get_voxel_doses.py ProductionPatients/WL/

#echo patient WT
#python lymphkill/structure_loading.py ProductionPatients/WT/
#python lymphkill/mask_generation.py ProductionPatients/WT/
#python get_voxel_doses.py ProductionPatients/WT/


#python pipeline.py

