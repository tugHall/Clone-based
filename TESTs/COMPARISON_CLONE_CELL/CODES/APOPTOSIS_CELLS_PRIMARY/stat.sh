#!/bin/sh
#QSUB -queue i9acc 
#QSUB -node 9     
#
#PBS -N bulkjob 
#PBS -l walltime=00:30:00 
#

date
cd ${PBS_O_WORKDIR}

###### mpijob ./a.out < ./param_3N.in > out

bulkjob   ./joblist.txt 

### R  --vanilla < SIM.R > out.job 

echo  
echo "##### done #####"

date

