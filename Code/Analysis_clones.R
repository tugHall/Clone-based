# Diversity of clones

library(stringr)

library(ape)
library(ggplot2)
library(ggtree)

par(xpd=TRUE, cex.lab=2, lwd = 2, mar = c(5, 5, 5, 5), tcl = 0.5, cex.axis = 1.75,  mgp = c(3, 0.6, 0))


source(file = "Code/tugHall_clone_functions.R")
onco <<- oncogene$new()        # make the vector onco about the hallmarks
onco$read(genefile)          # read the input info to the onco from genefile - 'gene_cds2.txt'
hall <<- hallmark$new()        # make a vector hall with hallmarks parameters
hall$read(genefile, onco$name)     # read from the genefile - 'gene_cds2.txt'

source("Code/Functions_clones.R")
# The function is in the Analyze.R file 
analyze_data(cloneoutfile)


ln <- strsplit(names(data_last), " ")
n_cluster <- which(substr(ln, 1, 9) == "PosDriver" )   # the numbers of columns for PosDriver to make clusters

# time_max - number of time steps 

######### Save CLONES: ##################
clones <- matrix(nrow = length(data_flow$Time), ncol = 6)
clones[,1] <- as.integer(as.character(data_flow$Time))
clones[,2] <- as.integer(as.character(data_flow$ID)) ## as.integer(as.character(data_flow$Clone.number))

x <- as.character(data_flow$ParentID.Birthday) 
x <- unlist(strsplit(x, split = ":") )
x <- as.integer(x)
clones[,3] <- x[c(TRUE,FALSE)]
clones[,4] <- x[c(FALSE,TRUE)]
clones[,5] <- as.integer(as.character(data_flow$N_cells))

dimnames(clones) <- list(1:length(data_flow$Time),c("Time","Clone_ID","Parent_ID","Birthday","Num_cells", "Cluster_ID") )

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
### the calculation of cluter's ID

sb_data <- subset.data.frame(data_flow,select = c(1, 4, n_cluster) )

sb_data[,1] <- as.integer(as.character(sb_data[,1]))
sb_data[,2] <- as.integer(as.character(sb_data[,2]))

for (i in 3:(2+onco$len)) {
  sb_data[,i] <- as.integer(sb_data[,i])
  sb_data[which(sb_data[,i] == 1),i] <- 0 
  sb_data[which(sb_data[,i] != 0),i] <- 1 
  
}

BinToDec <- function(x)  sum(2^(which(rev(unlist(strsplit(as.character(x), "")) == 1))-1))

sb_data[,(3+onco$len)] <- 0
names( sb_data ) <- c (names( sb_data[,1:(2+onco$len)] ), "Cluster_ID")

for (i in 1:length(sb_data[,1])) {

  inp <- paste0(sb_data[i,3:(2+onco$len)], collapse = "")
  sb_data[i,(3+onco$len)] <- BinToDec(inp)
  
  }


clones[,6] <- sb_data[,(3+onco$len)]
##################






diversity <- matrix(nrow = time_max, ncol = 2)

for (i in 1:time_max) {
diversity[i,1] <- length(clones[which(clones[,1]==i),2])   # number of clones
# diversity[i,2] <- diversity[i,1] / (data_avg$N[i] + data_avg$M[i])
 diversity[i,2] <- length(table(clones[which(clones[,1]==i),6])  )   # number of clusters
}  

# postscript(file =  "Figures/N_clones.eps", horizontal = FALSE, onefile = FALSE, paper = "special", height = 10, width = 10)

g_range_y <- range(1, max(diversity[,1]))

# plot(1:time_max,diversity[,1],type = "l",col = "red", lwd=2, ylim=g_range_y, xlab = "Time step", ylab = "Number of clones") # number of clones - diversity

plot(1:time_max,diversity[,1],type = "l", log = "y", col = "blue", lwd=2, ylim=g_range_y, xlab = "Time step", ylab = "Number of clones and clusters") ; axis(3, labels = FALSE) ; axis(4, labels = FALSE)

 lines(1:time_max,diversity[,2],type = "l",col = "red", lwd=2)

 legend(1, 1.15*g_range_y[2], c("Clones","Clusters"), 
       lwd=2,cex=1.5,col=c("blue","red"), lty = 1:1, horiz = TRUE, bty = "n")


save_fig("Figures/N_clones_clusters.eps")

rl <-  readline(prompt="This is a plot for Number of clones - Press Enter  ")

rm(sb_data)


x <- table(clones[,"Cluster_ID"]) # names of all clusters - "numeric" name of cluster

total_clusters <- length(names(x)) # number of all clones

 evolution_clusters <- data.frame(matrix(nrow = time_max, ncol = total_clusters))

 names(evolution_clusters) <- names(x)
 
for (i in 1:time_max) {
y <-  table(clones[which(clones[,1]==i),"Cluster_ID"])
evolution_clusters[i,match(names(y),names(evolution_clusters))] <- y  
}

 evolution_clusters[is.na(evolution_clusters)] <- 0


 tree_cluster <- calc_tree(evolution_clusters, clones, total_clusters)

 tree_cluster$label <- names(table(cluster_tree$Parent_ID))
 
 # ggtree(tree_cl) + theme_tree2()
 
 # label_nodes = names(table(clone_tree$Parent_ID))
 
 #p <- ggtree(tree_cl, color="blue", size=1.5, linetype=1)     + 
#   geom_nodepoint(mapping = NULL, data = NULL, position = "identity", na.rm = FALSE, show.legend = FALSE, color = "red", fill = "red", size = 3, shape=23)   +   
 #        geom_tippoint(color = "skyblue", size = 3) + 
  #         geom_tiplab(size = 6) + 
   #          geom_rootedge(color="blue", size=1.5, linetype=1) 

 
 
 p <- ggtree(tree_cluster, color="blue", size=1.5, linetype=1, ladderize = TRUE)     + 
   geom_nodepoint(mapping = NULL, data = NULL, position = "identity", na.rm = FALSE, show.legend = FALSE, color = "red", fill = "red", size = 3, shape=23)   +   
   geom_tippoint(color = "skyblue", size = 3)  
#    geom_rootedge(color="blue", size=1.5, linetype=1) 
 
 d <- p$data
 d$label[is.na(d$label)] <- tree_cluster$label
 
 
 plot( p + geom_text2(data=d, aes(label=label), nudge_x = 0.04, size = 6 )  )
 

 save_fig("Figures/ggtree_clusters.eps")
 
 rl <-  readline(prompt="This is a ggtree plot for  clusters - Press Enter  ") 
  ######## plot TREE:
  # plot(tree_cl, cex = 1.4)
  
 
  ### APE library:
  plotTreeTime(tree_cluster, tip.dates = cluster_tree$Time_start, show.tip.label = TRUE, label.offset = 0.01 )
  
  nodelabels(text = tree_cluster$label, frame = "circle")
  # tiplabels()
  # dgelabels()
  
  save_fig("Figures/Tree_clusters.eps")
 
  rl <-  readline(prompt="This is a Time tree plot for  clusters - Press Enter  , next can take a time, please, wait for a while - time depends on the number of all clones") 
  
  #dev.copy(pdf, "Figures/Tree_clones.pdf")    
  #dev.off()
  
  
  
  
  
  
  ########################################### CLONES: ##################################
  x <- table(clones[,"Clone_ID"]) # names of all clones - "numeric" name of clone
  
  total_clones <- length(names(x)) # number of all clones
  
  tree_clones <- calc_tree_clones(clones, total_clones)
  
  tree_clones$label <- names(table(clone_tree$Parent_ID))

  
  p <- ggtree(tree_clones, color="blue", size=0.5, linetype=1, ladderize = TRUE)  #   + 
  #  geom_nodepoint(mapping = NULL, data = NULL, position = "identity", na.rm = FALSE, show.legend = FALSE, color = "red", fill = "red", size = 3, shape=23)   +   
  #  geom_tippoint(color = "skyblue", size = 3)  
  #    geom_rootedge(color="blue", size=1.5, linetype=1) 
  
  #d <- p$data
  #d$label[is.na(d$label)] <- tree_clones$label
  
  plot( p ) # + geom_text2(data=d, aes(label=label), nudge_x = 0.04, size = 6 )  )
  
  save_fig("Figures/ggtree_clones.eps")
  
  rl <-  readline(prompt="This is a ggtree plot for clones - Press Enter  ") 
  ######## plot TREE:
  # plot(tree_cl, cex = 1.4)
  
  
  ### APE library:
  plotTreeTime(tree_clones, tip.dates = clone_tree$Time_start, show.tip.label = FALSE, label.offset = 0.01 )
  
  #nodelabels(text = tree_clones$label, frame = "circle")
  # tiplabels()
  # dgelabels()
  
  save_fig("Figures/Tree_clones.eps")
  
  rl <-  readline(prompt="This is a Time tree plot for  clones - Press Enter  ") 
  
  
  
  ########################################## END CLONES ######################################
  
  
  
  
  
  
  
  
  
  
  
  
  
# Make a large number of colors

nm <- length(evolution_clusters)

w <- (nm^(1/3)) %/% 1 +1

st <-  w^3 %/% nm

sq <- seq(0,1-1/w,1/w)

cr <- 1:nm

l <- 0
R <- 1:(w^3)
G <- R
B <- R

for (i in 1:w) {
  for (j in 1:w) {
    for (k in 1:w) {
      l <- l+1
      R[l] <- sq[i]
      G[l] <- sq[j]
      B[l] <- sq[k]
      
    } 
  }  
}

# seq(1,w^3,st) # is consequence of each color to make a high diversity of colors
jColor <- data.frame(number = 1:length(seq(1,w^3,st)),color = rgb(R[seq(1,w^3,st)],G[seq(1,w^3,st)],B[seq(1,w^3,st)]))

print("This is a plot for  number of cells in each cluster: ") 

for (ll in 1:2) {
  
  if (ll == 1) g_range_y <- range(1, max(evolution_clusters)) else g_range_y <- range(1,100)

  plot(1:time_max,evolution_clusters[,1],type = "l", log = "y", col = jColor$color[1],lwd=3,xlab = "Time step",
     ylab = "Number of cells in each cluster",ylim=g_range_y)

  for (i in 2:total_clusters) {
    lines(1:time_max,evolution_clusters[,i],col = jColor$color[i],lwd=3)
    }
  if (ll == 1) rl <-  readline("Press Enter for another scale")
 
  #if (ll == 1) dev.copy(pdf, "Figures/N_cells_in_clones_1.pdf")   else dev.copy(pdf, "Figures/N_cells_in_clones_2.pdf")    
  #dev.off()
  if (ll == 1)   save_fig("Figures/N_cells_in_clusters_1.eps")    else save_fig("Figures/N_cells_in_clusters_2.eps") 
  }

rl <-  readline("Press Enter")

# the clones at the last time step
cl <- clones[which(clones[,1]==time_max),]

# Most popular clone
barplot(table(cl[,"Cluster_ID"]),xlab = "The ID of cluster", ylab = "Number of cells in the cluster", log = "y",
         cex.name = 2, space=0.7, col = "green")

save_fig("Figures/Barplot_N_cells_in_clusters.eps") 
#dev.copy(pdf, "Figures/Barplot_N_cells_in_clones.pdf")    
#dev.off()
rl <-  readline(prompt="This is a barplot for Numbers of cells in clusters at last timestep - Press Enter ")


# Most popular clone (Drivers and Passengers)
barplot(sort(table(cl[,"Num_cells"]), decreasing = TRUE),xlab = "The ID of clone", ylab = "Number of cells in the clone",   log = "y", names.arg = "", #  main = "FOR DRIVERS AND PASSENGERS",
        cex.name = 1.6, space=0.7, col = "green")
save_fig("Figures/Barplot_N_cells_in_clones.eps")
#dev.copy(pdf, "Figures/Barplot_N_cells_in_clones_DP.pdf")    
#dev.off()
rl <-  readline(prompt="This is same barplot for DRIVERS AND PASSENGERS - Press Enter ")


########################################################################################################
###################################### INEQUALITY    COEFFICIENT  ######################################

ineq_clones <- matrix(0,nrow = (time_max+1),ncol = 2)

# Inequality measure 
for (k in 1:(time_max+1) ) {
  cl <- clones[which(clones[,"Time"]== (k-1) ),"Clone_ID"]
  ineq_clones[k,1] <- ineq(cl,type = "Gini")

  cl2 <- clones[which(clones[,"Time"]== (k-1) ),"Cluster_ID"]
  ineq_clones[k,2] <- ineq(cl2,type = "Gini")
}

ineq_clones[which(ineq_clones[,2] == "NaN"),] <-0 


g_range_y <- range(0, 1) # max(ineq_clones[,2], ineq_clones[,1]))

plot(0:time_max,ineq_clones[,1],type = "l", lwd=3, pch = 19, col = "blue", 
      ylim = g_range_y, xlab = "Time step", ylab = "Inequlity coefficient")
lines(0:time_max,ineq_clones[,2],type = "l", lwd=3, pch = 19, col =  "red")


# g_range_x[1]/2+g_range_x[2]/2.5, 1.2*g_range_y[2]
legend(1, 1.15*g_range_y[2], c("Clones","Clusters"), 
       lwd=2,cex=1.4,col=c("blue","red"), lty = 1:1, horiz = TRUE, bty = "n")

save_fig("Figures/Inequality.eps")
#dev.copy(pdf, "Figures/Inequality.pdf")    
#dev.off()
rl <-  readline(prompt="This is a plot for inequality coefficient - Press Enter ")



plot(data_avg$N,ineq_clones[,1],type = "l", lwd=3, pch = 19, col = "blue", ylim = g_range_y, xlab = "Number of primary tumor cells", ylab = "Inequlity coefficient")
points(data_avg$N,ineq_clones[,2],type = "l", lwd=3, pch = 19, col =  "red")
legend(min(data_avg$N)-1, 1.15*g_range_y[2], c("Clones","Clusters"),
       lwd=2,cex=1.4,col=c("blue","red"), lty = 1:1, horiz = TRUE, bty = "n")

save_fig("Figures/Inequality_primary.eps")
#dev.copy(pdf, "Figures/Inequality_primary.pdf")    
#dev.off()

rl <-  readline(prompt="This is a plot for inequality coefficient for primary tumor cells - Press Enter ")



plot(data_avg$M,ineq_clones[,1],type = "l", lwd=3, pch = 19, col = "blue", ylim = g_range_y, xlab = "Number of metastasis cells", ylab = "Inequlity coefficient")
points(data_avg$M,ineq_clones[,2],type = "l", lwd=3, pch = 19, col =  "red")
legend(min(data_avg$M), 1.15*g_range_y[2], c("Clones","Clusters"), 
       lwd=2,cex=1.4,col=c("blue","red"), lty = 1:1, horiz = TRUE, bty = "n")

save_fig("Figures/Inequality_metastasis.eps")
#dev.copy(pdf, "Figures/Inequality_metastasis.pdf")    
#dev.off()
rl <-  readline(prompt="This is a plot for inequality coefficient for metastasis cells - Press Enter ")



plot(data_avg$M+data_avg$N,ineq_clones[,1],type = "l", lwd=3, pch = 19, col = "blue", ylim = g_range_y, xlab = "Number of all cells", ylab = "Inequlity coefficient")
points(data_avg$M+data_avg$N,ineq_clones[,2],type = "l", lwd=3, pch = 19, col =  "red")
legend(min(data_avg$M+data_avg$N), 1.15*g_range_y[2], c("Clones","Clusters"),
       lwd=2,cex=1.4,col=c("blue","red"), lty = 1:1, horiz = TRUE, bty = "n")
save_fig("Figures/Inequality_all_cells.eps")
#dev.copy(pdf, "Figures/Inequality_all_cells.pdf")    
#dev.off()

rl <-  readline(prompt="This is a plot for inequality coefficient for all cells - Press Enter ")






################################## VAF 

### data_last has information for all mutations in genes of all cells 
### from 22 column         to 21 + length(onco) - for drivers and 
### from 22 + length(onco) to 21 + 2*length(onco) for passangers


VAF <- NULL
N <- data_last$N[1]
M <- data_last$M[1]
N_all <- N + M

st <- min(n_cluster)
for (k in st:(st-1 + 2 * onco$len)) {
  
  if (k > st-1+onco$len) DriverPasngr <- "P" else DriverPasngr <- "D" 
  if (k > st-1+onco$len) Gene <- onco$name[ k - st+1 - onco$len] else Gene <- onco$name[ k - st+1]
 
  d <- as.vector( data_last[,k] )
  d <- as.character(d)
  d <- str_replace_all(d,":.....,",",")
  d <- str_replace_all(d,":....,",",")
  d <- str_replace_all(d,":...,",",")
  d <- str_replace_all(d,":..,",",")
  d <- str_replace_all(d,":.,",",")

  d <- str_replace_all(d,":.....","")
  d <- str_replace_all(d,":....","")
  d <- str_replace_all(d,":...","")
  d <- str_replace_all(d,":..","")
  d <- str_replace_all(d,":.","")

  d_Primary    <- d[which(data_last$im  < 1)]
  d_Metastatic <- d[which(data_last$im == 1)]
  
  d <- str_split(d,",")  
  d <- unlist(d)
  d <- as.integer(d)

  d_Primary <- str_split(d_Primary,",")  
  d_Primary <- unlist(d_Primary)
  d_Primary <- as.integer(d_Primary)

  d_Metastatic <- str_split(d_Metastatic,",")  
  d_Metastatic <- unlist(d_Metastatic)
  d_Metastatic <- as.integer(d_Metastatic)  
  
  Z_M <- TRUE
  Z_N <- TRUE
  Z_M <- all(is.na(d_Metastatic))
  Z_N <- all(is.na(d_Primary))
  
  out <- NULL
  out <- as.data.frame( table(d) )
  if (!Z_M | !Z_N) names(out) <- c("pos","Freq")
  
  if (!Z_N && N > 0) {
      out_Prim <- NULL
      out_Prim <- as.data.frame( table(d_Primary) )
      names(out_Prim) <- c("pos","Freq_Prim")
      out <- merge.data.frame(out_Prim, out, by = "pos" , all = TRUE)
      
      } else out["Freq_Prim"] <- 0
  
  
  if (!Z_M && M > 0) {
      out_Met <- NULL
      out_Met <- as.data.frame( table(d_Metastatic) )
      names(out_Met) <- c("pos","Freq_Met")
      out <- merge.data.frame(out_Met, out, by = "pos" , all = TRUE)
      
      } else out["Freq_Met"] <- 0
  
  
      out[is.na(out)] <- 0 
  
      if ( N > 0 ) out$VAF_Prim <- 0.5 * out$Freq_Prim / N   else out$VAF_Prim <- 0
      if ( M > 0 ) out$VAF_Met  <- 0.5 * out$Freq_Met  / M   else out$VAF_Met  <- 0
      if (!Z_M | !Z_N) out$VAF <- 0.5 * out$Freq / N_all     else out$VAF      <- 0
      
            
  # nm <- names(data_last[k])
  VAF1 <- NULL
  VAF1 <- cbind.data.frame(DriverPasngr, Gene, out$pos, 
                out$VAF_Prim,  out$Freq_Prim, N, 
                out$VAF_Met ,  out$Freq_Met,  M,
                out$VAF,  out$Freq, N_all)
  VAF <- rbind.data.frame(VAF,VAF1)
  
  }

header <- c( "DriverPasngr", "Gene", "Position", 
             "VAF_Primary", "Ncells_Primary_wMutation", "Ncells_Primary",
             "VAF_Metastatic", "Ncells_Metastatic_wMutation", "Ncells_Metastatic", 
             "VAF_PriMet", "Ncells_PriMet_wMutation", "Ncells_PriMet" )

names(VAF) <- header
write.table(VAF,file = "Output/VAF.txt", append = FALSE, row.names = FALSE, sep="\t")

print("VAF is saved to the file `Output/VAF.txt` ")
