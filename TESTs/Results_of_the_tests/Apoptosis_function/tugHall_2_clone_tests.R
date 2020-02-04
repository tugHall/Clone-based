
##########################################################################################
###  The simulation uses the functions and classes in the "Code/tugHall_functions.R" 

library(stringr)

source(file = "Code/tugHall_clone_functions.R")


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

E0 <<- 1E-3       # parameter in the division probability  
F0 <<- 10        # parameter in the division probability  
m0 <<-  1E-9       # mutation probability  
uo <<- 0.5        # oncogene mutation probability  
us <<- 0.5        # suppressor mutation probability  
s0 <<-  10         # parameter in the sigmoid function  
k0 <<-  0        # Environmental death probability  

### Additional parameters of simulation
censore_n <<- 30000       # Max cell number where the program forcibly stops
censore_t <<- 2         # Max time where the program forcibly stops


##########################################################################################
### Simulation of the cancer cell/clone evolution:
model(genefile, clonefile, geneoutfile, cloneoutfile, logoutfile, E0, F0, m0, uo, us, s0, k0, censore_n, censore_t)



##########################################################################################
#### Analysis of the output data:

# See the output file: clonout.txt 


