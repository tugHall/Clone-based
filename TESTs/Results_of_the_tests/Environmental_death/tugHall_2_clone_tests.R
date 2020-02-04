
##########################################################################################
###  The simulation uses the functions and classes in the "Code/tugHall_functions.R" 

library(stringr)

source(file = "Code/tugHall_clone_functions_tests.R")


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

E0 <<- 1E-10       # parameter in the division probability  
F0 <<- 10        # parameter in the division probability  
m0 <<-  1E-9       # mutation probability  
uo <<- 0.5        # oncogene mutation probability  
us <<- 0.5        # suppressor mutation probability  
s0 <<-  10         # parameter in the sigmoid function  
k0 <<-  0.5        # Environmental death probability  

### Additional parameters of simulation
censore_n <<- 30000       # Max cell number where the program forcibly stops
censore_t <<- 100         # Max time where the program forcibly stops


##########################################################################################
# if you have a new format of gene file, please, use change of columns function like: 
# genefile <- changeCol(genefile)

### Making of the input file for initial clones

x <- 1
# For primary tumor cells: 
xz <- data.frame(V1=x,V2=as.character(""), V3=rep.int(1000,length(x)))   # the id of clone, the genes, the number of cells in the clone
#xz[1,2] <- "PIK3CA"
#xz[2,2] <- "APC"
#xz[3,2] <- "KRAS"
#xz[4,2] <- "TP53"
xz$V2 <- as.character(xz$V2)

write.table(xz,file = clonefile, col.names = FALSE,sep = "\t",row.names = FALSE)

##########################################################################################
### Simulation of the cancer cell/clone evolution:
# model(genefile, clonefile, geneoutfile, cloneoutfile, logoutfile, E0, F0, m0, uo, us, s0, k0, censore_n, censore_t)


# To repeat the 100 simulations and to mesure the time difference

Sys.setenv(TZ = "Japan")
print(Sys.timezone())

RP_AGAIN <- function()    {
  
  st <- Sys.time()
  
  print("The time at START:")
  print( Sys.time() )
  
  for (rp in 1:100 )  { 
    print(rp)
    model_test(genefile, clonefile, geneoutfile, cloneoutfile, logoutfile, E0, F0, m0, uo, us, s0, k0, censore_n, censore_t,rp)
  }
  
  print("The time at FINISH:")
  print( Sys.time() )
  nd <- Sys.time()
  print(nd - st)
  
}


RP_AGAIN()


##########################################################################################
#### Analysis of the output data:


source(file = "Code/Average_tests.R")
ave_plot(data_out = data_out, output = "N", lg = "Number of cells", xl = c(0,70), yl = c(0,1800), nm1 = "Fig_1.jpg", nm2 = "Fig_2.jpg")



####### SIMULATION FOR METASTASIS CELLS:

x <- 1
# For metastasis cells:
xz <- data.frame(V1=x,V2=as.character("PIK3CA,APC,KRAS,TP53"), V3=rep.int(1000,length(x)))   # the id of clone, the genes, the number of cells in the clone
xz$V2 <- as.character(xz$V2)

write.table(xz,file = clonefile, col.names = FALSE,sep = "\t",row.names = FALSE)


RP_AGAIN()


##########################################################################################
#### Analysis of the output data:
source(file = "Code/Average_tests.R")
ave_plot(data_out = data_out, output = "M", lg = "Number of cells", xl = c(0,105), yl = c(0,2000), nm1 = "Fig_3.jpg", nm2 = "Fig_4.jpg")










