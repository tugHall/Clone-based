
##########################################################################################
###  The simulation uses the functions and classes in the "Code/tugHall_2.1_functions.R" 

library(stringr)   # to use string data in input files
library(actuar)    # to use BIG NUMBERS in N_cell variable

source(file = "Code/tugHall_2.1_functions.R")


## Create folders:  /Input, /Output and /Figures 

mainDir <- getwd()
subDir <- "Output"
if (! file.exists(subDir)){  dir.create(file.path(mainDir, subDir)) }

subDir <- "Input"
if (! file.exists(subDir)){  dir.create(file.path(mainDir, subDir)) }

subDir <- "Figures"
if (! file.exists(subDir)){  dir.create(file.path(mainDir, subDir)) }


##########################################################################################
### Files to output and input data

genefile <- 'Input/gene_cds2.txt'    # gene file 
clonefile <- 'Input/cloneinit.txt'     # initial Cells 

### Output files
geneoutfile <- 'Output/geneout.txt'  # Gene Out file with Hallmarks 
cloneoutfile <- 'Output/cloneout.txt'  # output information of simulation
logoutfile <-  'Output/log.txt'      # log file to save the input information of simulation - "log.txt"
### Output/Weights.txt               # file with gene weights for hallmarks


##########################################################################################
# Probabilities of processes

E0 <<-  1E-3       # parameter in the division probability  
F0 <<-  2          # parameter in the division probability  
m0 <<-  1E-9       # mutation probability  
uo <<-  0.5        # oncogene mutation probability  
us <<-  0.5        # suppressor mutation probability  
s0 <<-  10         # parameter in the sigmoid function  
k0 <<-  0.2        # Environmental death probability  
d0 <<-  0   # Initial probability to divide cells
### Additional parameters of simulation
censore_n <<- 10^4       # Max cell number where the program forcibly stops
censore_t <<- 2         # Max time where the program forcibly stops


##########################################################################################
# if you have a new format of gene file, please, use change of columns function like: 
# genefile <- changeCol(genefile)

##########################################################################################
### Simulation of the cancer cell/clone evolution:
model(genefile, clonefile, geneoutfile, cloneoutfile, logoutfile, E0, F0, m0, uo, us, s0, k0, censore_n, censore_t, d0)

##########################################################################################
#### Analysis of the output data:


### Please, just to see cloneout file in the Output folder and compare data with template.

