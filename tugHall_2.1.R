
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

E0 <<-  1E-4       # parameter in the division probability  
F0 <<-  10         # parameter in the division probability  
m0 <<-  1E-7       # mutation probability  
uo <<-  0.5        # oncogene mutation probability  
us <<-  0.5        # suppressor mutation probability  
s0 <<-  10         # parameter in the sigmoid function  
k0 <<-  0.2        # Environmental death probability  
d0 <<-  0.35   # Initial probability to divide cells
### Additional parameters of simulation
censore_n <<- 10^5       # Max cell number where the program forcibly stops
censore_t <<- 100         # Max time where the program forcibly stops


##########################################################################################
# if you have a new format of gene file, please, use change of columns function like: 
# genefile <- changeCol(genefile)

### Making of the input file for initial clones

x <- 1
#xz <- data.frame(V1=x,V2=as.character("PIK3CA,APC,KRAS,TP53"), V3=rep.int(1000,length(x)))   # the id of clone, the genes, the number of cells in the clone
xz <- data.frame(V1=x,V2=as.character(""), V3=rep.int(10^3,length(x)))   # the id of clone, the genes, the number of cells in the clone
#xz[1,2] <- "PIK3CA"
xz[2,] <- c(2,"APC",10^3)
#xz[3,2] <- "KRAS"
#xz[4,2] <- "TP53"
xz$V2 <- as.character(xz$V2)

write.table(xz,file = clonefile, col.names = FALSE,sep = "\t",row.names = FALSE)

##########################################################################################
### Simulation of the cancer cell/clone evolution:
model(genefile, clonefile, geneoutfile, cloneoutfile, logoutfile, E0, F0, m0, uo, us, s0, k0, censore_n, censore_t, d0)

##########################################################################################
#### Analysis of the output data:







# Note: if output files have no data to plot Code/Analysis.R produces errors during plotting

source("Code/Analysis_clones.R")



##########################################################################################

# In order to make report, please, use USER-GUIDE.Rmd to show results of simulation

# In order to improve the output plot, please, use Code/Functions.R and Code/Analysis.R scripts

