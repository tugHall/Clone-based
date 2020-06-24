############################# Libraries ######################################

library(stringr)
library(ape)
library(ggplot2)
library(ggtree)

par(xpd=TRUE, cex.lab=2, lwd = 2, mar = c(5, 5, 5, 5), tcl = 0.5, cex.axis = 1.75,  mgp = c(3, 0.6, 0))

#################### Diversity of clones  ####################################
source(file = "Code/tugHall_2.1_functions.R")

onco <<- oncogene$new()        # make the vector onco about the hallmarks
onco$read(genefile)          # read the input info to the onco from genefile - 'gene_cds2.txt'
hall <<- hallmark$new()        # make a vector hall with hallmarks parameters
hall$read(genefile, onco$name)     # read from the genefile - 'gene_cds2.txt'

source("Code/Functions_clones.R")
# The function is in the Code/Functions_clones.R file 
analyze_data(cloneoutfile)

ln <- strsplit(names(data_last), " ")
n_cluster <- which(substr(ln, 1, 9) == "PosDriver" )   # the numbers of columns for PosDriver to make clusters

# time_max - number of time steps 

clones <- Make_clones(data_flow = data_flow)
##################

Pl_diversity(time_max = time_max, clones = clones)

Pl_trees(time_max = time_max, clones = clones)

########################################## END CLONES ######################################
  
  
########################################## CELLS IN CLONES and CLUSTERs ###############################
#################################### Number of cells in CLUSTERs   ####################################
 
evolution_clusters <- Pl_number_cells_in_clusters(clones = clones)

evolution_clones  <-  Pl_number_cells_in_clones(clones = clones)

# the clones at the last time step
cl <- evolution_clones[time_max,]
cl <- cl[which(cl != 0)]
cl.df <- data.frame(Variable = names(cl), Count = as.numeric(t(cl)) )

cs <- evolution_clusters[time_max,]
cs <- cs[which(cs != 0)]
cs.df <- data.frame(Variable = names(cs), Count = as.numeric(t(cs)) )


b_ggpl(data.sim = cl.df, x_l = "ID of clones", y_l = "Number of cells")
save_fig("Figures/Barplot_N_cells_in_clones.eps") 
rl <-  readline(prompt="This is a barplot for Numbers of cells in clones at last timestep - Press Enter ")

b_ggpl(data.sim = cs.df, x_l = "ID of cluster", y_l = "Number of cells")
save_fig("Figures/Barplot_N_cells_in_clusters.eps") 
rl <-  readline(prompt="This is a barplot for Numbers of cells in clusters at last timestep - Press Enter ")


########################################################################################################
###################################### INEQUALITY    COEFFICIENT  ######################################

ineq_clones <- ineq_time(time_max = time_max, clones = clones)

g_range_y <- range(0, max(ineq_clones) )
plot(0:time_max,ineq_clones,type = "l", lwd=3, pch = 19, col = "blue", 
      ylim = g_range_y, xlab = "Time step", ylab = "Inequlity coefficient")

save_fig("Figures/Inequality.eps")
rl <-  readline(prompt="This is a plot for inequality coefficient - Press Enter ")

plot(data_avg$N,ineq_clones,type = "l", lwd=3, pch = 19, col = "blue", ylim = g_range_y, xlab = "Number of primary tumor cells", ylab = "Inequlity coefficient")
#points(data_avg$N,ineq_clones[,2],type = "l", lwd=3, pch = 19, col =  "red")
#legend(min(data_avg$N)-1, 1.15*g_range_y[2], c("Clones","Clusters"),
#       lwd=2,cex=1.4,col=c("blue","red"), lty = 1:1, horiz = TRUE, bty = "n")

save_fig("Figures/Inequality_primary.eps")
rl <-  readline(prompt="This is a plot for inequality coefficient for primary tumor cells - Press Enter ")

plot(data_avg$M,ineq_clones,type = "l", lwd=3, pch = 19, col = "blue", ylim = g_range_y, xlab = "Number of metastasis cells", ylab = "Inequlity coefficient")
#points(data_avg$M,ineq_clones[,2],type = "l", lwd=3, pch = 19, col =  "red")
#legend(min(data_avg$M), 1.15*g_range_y[2], c("Clones","Clusters"), 
#       lwd=2,cex=1.4,col=c("blue","red"), lty = 1:1, horiz = TRUE, bty = "n")

save_fig("Figures/Inequality_metastasis.eps")
rl <-  readline(prompt="This is a plot for inequality coefficient for metastasis cells - Press Enter ")

plot(data_avg$M+data_avg$N,ineq_clones,type = "l", lwd=3, pch = 19, col = "blue", ylim = g_range_y, xlab = "Number of all cells", ylab = "Inequlity coefficient")
#points(data_avg$M+data_avg$N,ineq_clones[,2],type = "l", lwd=3, pch = 19, col =  "red")
#legend(min(data_avg$M+data_avg$N), 1.15*g_range_y[2], c("Clones","Clusters"),
#       lwd=2,cex=1.4,col=c("blue","red"), lty = 1:1, horiz = TRUE, bty = "n")

save_fig("Figures/Inequality_all_cells.eps")
rl <-  readline(prompt="This is a plot for inequality coefficient for all cells - Press Enter ")

################################## VAF ##################################

### data_last has information for all mutations in genes of all cells 
### from 23 column         to 22 + length(onco) - for drivers and 
### from 23 + length(onco) to 22 + 2*length(onco) for passangers
VAF <- Safe_VAF(data_last = data_last, n_cluster = n_cluster)

