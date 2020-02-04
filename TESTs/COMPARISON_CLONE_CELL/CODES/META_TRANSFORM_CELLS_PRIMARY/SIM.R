

##########################################################################################
###  The simulation uses the functions and classes in the "Code/tugHall_functions.R" 

source(file = "Code/tugHall_functions.R")
### source("Code/Functions.R")

library(stringr)    # to calculate VAF

## Create folders:  /Input, /Output 

mainDir <- getwd()
subDir <- "Output"
if (! file.exists(subDir)){  dir.create(file.path(mainDir, subDir)) }

subDir <- "Input"
if (! file.exists(subDir)){  dir.create(file.path(mainDir, subDir)) }


##########################################################################################
### Files to output and input data

genefile <- 'Input/gene_cds2.txt'    # gene file 
cellfile <- 'Input/cellinit.txt'     # initial Cells 

### Output files
geneoutfile <- 'geneout.txt'  # Gene Out file with Hallmarks 
celloutfile <- 'cellout.txt'  # output information of simulation
logoutfile <-  'log.txt'      # log file to save the input information of simulation - "log.txt"
### Output/Weights.txt               # file with gene weights for hallmarks

### Making of the input file for initial NORMAL cells
x <- 1:1000
xz <- data.frame(x,"")    # (x, "PIK3CA,APC,KRAS,TP53" )   ### for metstasis
write.table(xz,file = cellfile, col.names = FALSE,sep = "\t",row.names = FALSE)

############### FUNCTION FOR ITTERATION ########################

SIMULATION = function(i, mainDir, numCores, onco, hall, genefile, cellfile, geneoutfile, 
               celloutfile, logoutfile  )  {
    source(file = "Code/tugHall_functions.R")
    # source("Code/Functions.R")

    # library(stringr)    # to calculate VAF

    ## Create folders:  /Input, /Output 

    mainDir <- getwd()
    subDir <- "Output"
    if (! file.exists(subDir)){  dir.create(file.path(mainDir, subDir)) }

    subDir <- "Input"
    if (! file.exists(subDir)){  dir.create(file.path(mainDir, subDir)) }


    ##########################################################################################
    ### Files to output and input data

    genefile <- 'Input/gene_cds2.txt'    # gene file 
    cellfile <- 'Input/cellinit.txt'     # initial Cells 

    ### Output files
    geneoutfile <- 'geneout.txt'  # Gene Out file with Hallmarks 
    celloutfile <- 'cellout.txt'  # output information of simulation
    logoutfile <-  'log.txt'      # log file to save the input information of simulation - "log.txt"
    ### Output/Weights.txt               # file with gene weights for hallmarks

    ## Define initial halmarks and genes and weights (we will change it later)

    onco <<- oncogene$new()        # make the vector onco about the hallmarks
    onco$read(genefile)          # read the input info to the onco from genefile - 'gene_cds2.txt'
    hall <<- hallmark$new()        # make a vector hall with hallmarks parameters
    hall$read(genefile, onco$name)     # read from the genefile - 'gene_cds2.txt'
        
    subDir <- paste0("Output/",i)
    if (! file.exists(subDir)){  dir.create(file.path(mainDir, subDir)) }
        
    ##########################################################################################
    # Probabilities of processes
        
    E0 <<- 1E-10       # parameter in the division probability  
    m <<-  1E-9       # mutation probability  
    k <<-  0   ## sample(c(0.1, 0.2,0.3,0.4,0.5, 0.6,0.7,0.8, 0.9),1,replace = TRUE)          # Environmental death probability 
    F0 <<- 10  ## sample(c(1,10,100,1000),1,replace = TRUE)        # parameter in the division probability 
    s <<-  10  ## sample(c(10,20,50,90),1,replace = TRUE)             # parameter in the sigmoid function  
    uo <<- 0.5 ## sample(c(0.01,0.1,0.2,0.3,0.4,0.5),1,replace = TRUE)         # oncogene mutation probability  
    us <<- 0.5 ## sample(c(0.01,0.1,0.2,0.3,0.4,0.5),1,replace = TRUE)         # suppressor mutation probability  
 
    ### Additional parameters of simulation
    censore_n <<- 30000       # Max cell number where the program forcibly stops
    censore_t <<- 100         # Max time where the program forcibly stops
        
    ### Write input data
    write_geneout(paste0("Output/",i,"/",geneoutfile), hall)                   # write the geneout.txt file with initial hallmarks 
        
    ######################## SIMULATION ########################
        
    model(onco, hall, genefile, cellfile, geneoutfile,  paste0("Output/",i,"/",celloutfile), logoutfile, E0, F0, m, uo, us, s, k, censore_n, censore_t)
        
       
    file_Cell_Out <- paste0("Output/",i,"/",celloutfile)

    # if ( file.exists( file_Cell_Out ) ) calc_VAF(i, file_Cell_Out )
}


############### START OF ITTERATION ########################

 SIM <- function(i){
   SIMULATION(i, mainDir, numCores, onco, hall, genefile, cellfile, geneoutfile, celloutfile, logoutfile) 
 }

#################### Parallel implementation ###############

# library(doParallel)
 library(parallel)
# library(foreach)

 numCores <- detectCores()
 print(numCores)
 
 node <- 1 
 N_sim <- 120   # numCores * 25  # number of simulation in 1 node

 inputs = ( N_sim*(node - 1) + 1):( N_sim*(node - 1) + N_sim )
 print(inputs)
 
 mclapply(inputs, SIM, mc.cores = numCores)

########################################################


