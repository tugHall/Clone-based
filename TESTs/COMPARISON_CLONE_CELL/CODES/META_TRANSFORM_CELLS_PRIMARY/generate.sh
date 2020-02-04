#!/bin/sh
#

date

###  R  --vanilla < SIM.R > out.job 

### let dance with 1000 simulations!!! 


### sed -i 's//new/g' file.txt

### sed "s/node <- 1/node <- $i/g" SIM.R > SIM_$i.R 

## remove from previus calculations
rm joblist.txt 


for ((i=1; i < 10; i++))
do

### make SIM_XXX.R file for R simultion 
sed "s/node <- 1/node <- $i/g" SIM.R > SIM_$i.R

### make the submission file stat_XXX.sh for HPC system

### sed "s/SIM.R > out.job/SIM_$i.R > out_$i.job/g" stat.sh > stat_$i.sh

echo   "R  --vanilla < SIM_$i.R > out_$i.job;1;S;  "  >> joblist.txt  


### qsub stat_$i.sh 

## f="PATHS/path_"$i".dat"

## make a simulation 

### cp path.dat $f

echo $i

done

echo  
echo "##### done #####"

date

