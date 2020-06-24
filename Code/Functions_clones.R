# The code to analize the results of simulation

save_fig <- function(file_input) {
    file_input <- paste0(substr(file_input,1,nchar(file_input)-4),".jpg")
    #dev.copy2eps(file = file_input, height = 10, width = 10) 
    dev.copy(jpeg,file_input, width = 10, height = 10, units = 'in', res = 300)
    dev.off()
    par(xpd=TRUE, cex.lab=2, lwd = 2, mar = c(5, 5, 5, 5), tcl = 0.5, cex.axis = 1.75,  mgp = c(3, 0.6, 0))
}


jCol <- function(nm = 12){
    
    # Make a large number of colors
    
    # nm <- length(evolution_clones)
    
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
    
    return(jColor)
  
}


b_ggpl <- function(data.sim, x_l = "x axis", y_l = "y axis") {
    
    g <- ggplot(data = data.sim, aes(x = Variable, y = Count), xlab = "" )
    qb <- geom_bar(stat = "identity", width=0.7)
    
    ax <- theme(axis.text.x = element_text(face="bold",      # color="#993333", 
                                           size=14, angle=0),
                axis.text.y = element_text(face="bold",      # color="#993333", 
                                           size=14, angle=0),
                axis.title.x = element_text(face="bold",       color="#FF3333", 
                                            size=24, angle=0),
                axis.title.y = element_text(face="bold",       color="#FF3333", 
                                            size=24, angle=90),
    )
    ax1 <- scale_x_discrete(x_l, labels=data.sim$Variable) 
    ax2 <- scale_y_continuous(name = y_l) 
    
    print(g + qb + ax + ax1 + ax2 )
}


######### Save CLONES: ##################
BinToDec <- function(x)  sum(2^(which(rev(unlist(strsplit(as.character(x), "")) == 1))-1))

Make_clones <- function( data_flow = data_flow ) {
  
    clones <- matrix(nrow = length(data_flow$Time), ncol = 6)
    clones[,1] <- as.integer(as.character(data_flow$Time))
    clones[,2] <- as.integer(as.character(data_flow$ID)) ## as.integer(as.character(data_flow$Clone.number))
    
    x <- as.character(data_flow$ParentID.Birthday) 
    x <- unlist(strsplit(x, split = ":") )
    x <- as.integer(x)
    clones[,3] <- x[c(TRUE,FALSE)]
    clones[,4] <- x[c(FALSE,TRUE)]
    clones[,5] <- as.numeric(as.character(data_flow$N_cells))
    
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
      sb_data[is.na(sb_data[,i]),i] <- 0     
    }
    
    sb_data[,(3+onco$len)] <- 0
    names( sb_data ) <- c (names( sb_data[,1:(2+onco$len)] ), "Cluster_ID")
    
    for (i in 1:length(sb_data[,1])) {
      inp <- paste0(sb_data[i,3:(2+onco$len)], collapse = "")
      sb_data[i,(3+onco$len)] <- BinToDec(inp)
    }
    
    clones[,6] <- sb_data[,(3+onco$len)]
    
    return(clones)
}



Pl_diversity <- function(time_max = time_max, clones = clones) {
    
    diversity <- matrix(nrow = time_max, ncol = 2)
    
    for (i in 1:time_max) {
      diversity[i,1] <- length(clones[which(clones[,1]==i),2])   # number of clones
      # diversity[i,2] <- diversity[i,1] / (data_avg$N[i] + data_avg$M[i])
      diversity[i,2] <- length(table(clones[which(clones[,1]==i),6])  )   # number of clusters
    }  
    
    g_range_y <- range(1, max(diversity[,1]))
    
    plot(1:time_max,diversity[,1],type = "l", log = "y", col = "blue", lwd=2, ylim=g_range_y, xlab = "Time step", ylab = "Number of clones and clusters") ; axis(3, labels = FALSE) ; axis(4, labels = FALSE)
    
    lines(1:time_max,diversity[,2],type = "l",col = "red", lwd=2)
    
    legend(1, 1.15*g_range_y[2], c("Clones","Clusters"), 
           lwd=2,cex=1.5,col=c("blue","red"), lty = 1:1, horiz = TRUE, bty = "n")
    
    save_fig("Figures/N_clones_clusters.eps")
    
    rl <-  readline(prompt="This is a plot for Number of clones - Press Enter  ")
}



Pl_trees <- function(time_max = time_max, clones = clones) {
  
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
    
    rl <-  readline(prompt="This is a Time tree plot for  clusters - Press Enter ") 
    
    
    ########################################### CLONES: ##################################
    x <- table(clones[,"Clone_ID"]) # names of all clones - "numeric" name of clone
    
    total_clones <- length(names(x)) # number of all clones
    
    if (total_clones > 100) {total_clones <- 100 }
    
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
    
    rl <-  readline(prompt="This is a ggtree plot for clones (or for first 100 clones) - Press Enter  ") 
    ######## plot TREE:
    # plot(tree_cl, cex = 1.4)
    
    if (  min(clone_tree$Time_start) == max(clone_tree$Time_start)  ) {  clone_tree$Time_start[2:total_clones] <- clone_tree$Time_start[2:total_clones] + 0.1 }
    ### APE library:
    plotTreeTime(tree_clones, tip.dates = clone_tree$Time_start, show.tip.label = FALSE, label.offset = 0.01 )
    
    
    #nodelabels(text = tree_clones$label, frame = "circle")
    # tiplabels()
    # dgelabels()
    
    save_fig("Figures/Tree_clones.eps")
    
    rl <-  readline(prompt="This is a Time tree plot for  clones (or for first 100 clones) - Press Enter  ") 
  
}  



Pl_number_cells_in_clusters <- function(clones = clones) {
  
    x <- table(clones[,"Cluster_ID"]) # names of all clusters - "numeric" name of cluster
    
    total_clusters <- length(names(x)) # number of all clones
    
    evolution_clusters <- data.frame(matrix(nrow = time_max, ncol = total_clusters))
    
    names(evolution_clusters) <- names(x)
    
    for (i in 1:time_max) {
        y <-  clones[which(clones[,1]==i),"Num_cells"]
        x <-  clones[which(clones[,1]==i),"Cluster_ID"]
        nm <- names( table(x) )
        for (k in 1:length(nm))    evolution_clusters[i, nm[k] ] <- sum(y[which( x == nm[k] )])
    }
    
    evolution_clusters[is.na(evolution_clusters)] <- 0
    
    jColor <- jCol( nm = length(evolution_clusters) )
    
    print("This is a plot for  number of cells in each cluster: ") 
    
    for (ll in 1:2) {
        if (ll == 1) {
            g_range_y <- range(1, max(evolution_clusters)) 
            sc <- "y"
        } else {
            g_range_y <- range(1,100)
            sc <- ""
      }
      
      plot(1:time_max,evolution_clusters[,1],type = "l", log = sc, col = jColor$color[1],lwd=3,xlab = "Time step",
           ylab = "Number of cells in each cluster",ylim=g_range_y)
      
      for (i in 2:total_clusters) {
         lines(1:time_max,evolution_clusters[,i],col = jColor$color[i],lwd=3)
      }
      
      if (ll == 1) { 
          rl <-  readline("Press Enter for another scale")
          save_fig("Figures/N_cells_in_clusters_1.eps")    
       }  else  {
          save_fig("Figures/N_cells_in_clusters_2.eps") 
      }
    }
    rl <-  readline("Press Enter")
    return(evolution_clusters)
}    




Pl_number_cells_in_clones <- function(clones = clones) {
      
    x <- table(clones[,"Clone_ID"]) # names of all clusters - "numeric" name of cluster
    
    total_clones <- length(names(x)) # number of all clones
    if (total_clones > 100)    total_clones <- 100 
    
    x <- x[1:total_clones]
    evolution_clones <- data.frame(matrix(data =0, nrow = time_max, ncol = total_clones))
    
    names(evolution_clones) <- names(x)
    
    for (i in 1:time_max) {
        y <-  clones[which(clones[,1]==i),"Num_cells"]
        x <-  clones[which(clones[,1]==i),"Clone_ID"]
        r <- which(x < total_clones + 1)
        x <- x[r]
        y <- y[r]
        nm <- names( table(x) )
        for (k in 1:length(nm))    evolution_clones[i, nm[k] ] <- sum(y[which( x == nm[k] )])
    }
    
    evolution_clones[is.na(evolution_clones)] <- 0
    
    jColor <- jCol( nm = length(evolution_clones) )
    
    print("This is a plot for  number of cells in each cluster: ") 
    
    for (ll in 1:2) {
        if (ll == 1) {
            g_range_y <- range(1, max(evolution_clones)) 
            sc <- "y"
        } else {
            g_range_y <- range(1,100)
            sc <- ""
      }
      
      plot(1:time_max,evolution_clones[,1],type = "l", log = sc, col = jColor$color[1],lwd=3,
           xlab = "Time step",
           ylab = "Number of cells in each clone",ylim=g_range_y)
      
      for (i in 2:total_clones) {
          lines(1:time_max,evolution_clones[,i],col = jColor$color[i],lwd=3)
      }
      
      if (ll == 1) { 
          rl <-  readline("Press Enter for another scale")
          save_fig("Figures/N_cells_in_clones_1.eps")    
      }  else  {
          save_fig("Figures/N_cells_in_clones_2.eps") 
      }
    }
    rl <-  readline("Press Enter")
    return(evolution_clones)
}    


ineq_time <- function(time_max,clones) {    
    ineq_clones <- matrix(0,nrow = (time_max+1),ncol = 1)
  
    # Inequality measure 
    for (k in 1:(time_max+1) ) {
        cl.x <- clones[which(clones[,"Time"]== (k-1) ),"Clone_ID"]
        cl.y <- clones[which(clones[,"Time"]== (k-1) ),"Num_cells"]
        if (max(cl.y) < 10^6) {
            cl <- rep(cl.x,cl.y)
            ineq_clones[k,1] <- ineq(cl,type = "Gini")
            } else {
                ineq_clones[k,1] <- 1
        }
    }
  
    return(ineq_clones)
}    


tab_clones <- function(d , freq){
    
    freq <- as.numeric(as.character(freq))
    tb <- table(d)
    for (h in 1:length(tb)) {
        pos <- as.integer(names(tb)[h])
        tb[h] <- sum(freq[which(d == pos)])
    }
    
    return(tb)
}



Safe_VAF <- function(data_last, n_cluster) {
    
    VAF <- NULL
    N <- data_last$N[1]
    M <- data_last$M[1]
    N_all <- N + M
    
    st <- min(n_cluster)
    print(paste0("Start column of genes is  ",st))
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
        
        # condition for Primary and metastasis cells:
        cond_prim <- which(data_last$im  < 1)
        cond_meta <- which(data_last$im == 1)
          
          
        d_Primary    <- d[cond_prim]
        d_Metastatic <- d[cond_meta]
        
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
        
        if (!Z_M | !Z_N) {
            out <- as.data.frame( tab_clones(d = d, freq = data_last$N_cells) )
            names(out) <- c("pos","Freq")
        } else {
            out <- data.frame(pos = 0, Freq = 0)
        }
        
        
        if (!Z_N && N > 0) {
            out_Prim <- NULL
            out_Prim <- as.data.frame( tab_clones(d = d_Primary, freq = data_last$N_cells[cond_prim]) )
            names(out_Prim) <- c("pos","Freq_Prim")
            out <- merge.data.frame(out_Prim, out, by = "pos" , all = TRUE)
          
        } else out["Freq_Prim"] <- 0
        
        
        if (!Z_M && M > 0) {
            out_Met <- NULL
            out_Met <- as.data.frame( tab_clones(d = d_Metastatic, freq = data_last$N_cells[cond_meta]) )
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
    
    ### Exclude this part, because of debugging for many simulations, 
    ### so VAF file can content the NA and zero data, if there is no mutations in some genes.
    
    # VAF <- VAF[-which(is.na.data.frame(VAF$Position)) , ]
    # row.names(VAF) <- 1:length(VAF$DriverPasngr)
    
    
    write.table(VAF,file = "Output/VAF.txt", append = FALSE, row.names = FALSE, sep="\t")
    
    print("VAF is saved to the file `Output/VAF.txt` ")
    
    return(VAF)
}





analyze_data <- function(cloneoutfile) {
    # Analising of results:
    data_out <- read.csv(cloneoutfile, sep="\t")
    data_out[is.na(data_out)] <- ""
    # make a readible names
    names(data_out)[6] <-  "c"
    names(data_out)[7] <-  "d"
    names(data_out)[8] <-  "i"
    names(data_out)[9] <-  "im"
    names(data_out)[10] <-  "a"
    names(data_out)[11] <- "k"
    names(data_out)[12] <- "E"
    names(data_out)[14] <- "Nmax"
    names(data_out)
  
    # average data
    data_avg <<- data_out[which(data_out$AvgOrIndx == "avg"),]
  
    # data without averaging - flow data
    data_flow <<- data_out[which(!data_out$AvgOrIndx == "avg"),]
    
    # the data of the last time step 
    time_max <<- max(data_flow$Time)
    data_last <<- data_flow[which(data_flow$Time == time_max),]
    rm(data_out)
    
          
            # let draw graphics 
            # Numbers of Metastasis and primary tumor cells
            g_range_y <- range(0, data_avg$N,data_flow$M)
            g_range_x <- range(min(data_avg$Time),max(data_flow$Time))
            plot(data_avg$Time,data_avg$M,type = "l",xlab = "Time step",
             ylab = "Number of cells",ylim=g_range_y,xlim = g_range_x,col = "red") ; axis(3, labels = FALSE) ; axis(4, labels = FALSE)
            lines(data_avg$Time,data_avg$N,type = "l",col = "blue")
            # g_range_x[1]/2+g_range_x[2]/2.5, 1.2*g_range_y[2]
            legend(g_range_x[1], 1.15*g_range_y[2], c("Primary tumor","Metastasis"), 
               lwd=2,cex=1.3,col=c("blue","red"), lty = 1:1,horiz = TRUE, bty = "n")
    #dev.copy(pdf, "Figures/N_cells.pdf")    
    #dev.off()
  
            save_fig("Figures/N_cells.eps")
    
    rl <-  readline(prompt="This is a plot for Numbers of primary tumor and Metastasis cells - Press Enter  ")
  
  
    # Average values of probabilities
    g_range_y <- range(0, data_avg[7:11])
    g_range_x <- range(min(data_avg$Time),max(data_avg$Time))
    
        
          plot(data_avg$Time,data_avg$d,type = "l", ylim=g_range_y,xlim = g_range_x,
           xlab = "Time step",ylab = "Average probabilities",col = "red") ; axis(3, labels = FALSE) ; axis(4, labels = FALSE)
      
          lines(data_avg$Time,data_avg$i,type = "l",col = "blue")
          lines(data_avg$Time,data_avg$im,type = "l",col = "green")
          lines(data_avg$Time,data_avg$a,type = "l",col = "orange")
          lines(data_avg$Time,data_avg$k,type = "l",col = "black")
  
          legend(g_range_x[1], 1.15*g_range_y[2], c("d","i","im","a","k"), cex=1.3,
                 col=c("red","blue","green","orange","black"), lty = 1:1,lwd = 2,horiz = TRUE, bty = "n")
        
  
          #dev.copy(pdf, "Figures/Probabilities.pdf")    
          #dev.off()
          
          save_fig("Figures/Probabilities.eps")
          
    rl <-  readline(prompt="This is a plot for Average values of probabilities - Press Enter  ")
  
      # The averaged values of Hallmarks 
      g_range_y <- range(0, data_avg[16:20])
      g_range_x <- range(min(data_avg$Time),max(data_avg$Time))
    
      plot(data_avg$Time,data_avg$Hd,type = "l", ylim=g_range_y,xlim = g_range_x,
         xlab = "Time step",ylab = "Averaged Hallmarks values",col = "red") ; axis(3, labels = FALSE) ; axis(4, labels = FALSE)
    
      lines(data_avg$Time,data_avg$Hi,type = "l",col = "blue",lwd = 2)
      lines(data_avg$Time,data_avg$Him,type = "l",col = "green",lwd = 2)
      lines(data_avg$Time,data_avg$Ha,type = "l",col = "orange",lwd = 2)
      lines(data_avg$Time,data_avg$Hb,type = "l",col = "black",lwd = 2)
  
      legend(g_range_x[1], 1.15*g_range_y[2], c("Hd","Hi","Him","Ha","Hb"), cex=1.2,
             col=c("red","blue","green","orange","black"), lty = 1:1,lwd = 2,horiz = TRUE, bty = "n")
  
  
      #dev.copy(pdf, "Figures/Hallmarks.pdf")    
      #dev.off()
      save_fig("Figures/Hallmarks.eps")  
    
    rl <-  readline(prompt="This is a plot for Average values of Hallmarks - Press Enter, but next takes a time for calculations")
  
    onco <<- oncogene$new()        # make the vector onco about the hallmarks
    onco$read(genefile)          # read the input info to the onco from genefile - 'gene_cds2.txt'
    
    order_dysfunction <<- data.frame(data_last[,22:(22+length(onco$name) + 1)],stringsAsFactors = FALSE)
    names(order_dysfunction) <<- c("ID",onco$name, "Num_cells")
    order_dysfunction[,1] <<- data_last$ID
    order_dysfunction[,length(onco$name) + 2] <<- data_last$N_cells
    
    for (i in 1:(length(onco$name)+2)) {
      order_dysfunction[,i] <<- as.character(order_dysfunction[,i])
    }
    #str(order_dysfunction)
    #order_dysfunction$ID <- as.character(order_dysfunction$ID)
    order_dysfunction$Num_cells <<- as.numeric(order_dysfunction$Num_cells)
    #str(order_dysfunction)
  
  
  
    # substr(x, regexpr(":",x)+1, ifelse(isTRUE(grep(",",x)>0),(regexpr(",",x)-1),nchar(x)))
  
  
    for (k in 1:length(order_dysfunction[,1])) {
  
        for (i in 2:(length(onco$name)+1)) {
  
            x<- order_dysfunction[k,i]
           order_dysfunction[k,i] <<- substr(x, regexpr(":",x)+1, ifelse(isTRUE(grep(",",x)>0),(regexpr(",",x)-1),nchar(x)))
        }
    }
  
  
    for (i in 1:length(onco$name)+1) {
        order_dysfunction[,i] <<- as.integer(order_dysfunction[,i])
    }
  
    outfile <- 'Output/Order_of_dysfunction.txt'
    header <- c('Order of gene dysfunction: from first to last', 'Number of clones with same order' , 'Number of cells with same order')
    write(header, outfile, append=FALSE, ncolumn=length(header), sep="\t")  
    
    #for (i in 1:length(order_dysfunction[,1])) {
    #  data <- c(order_dysfunction$ID[i],names(sort(order_dysfunction[i,2:(length(onco$name)+1)])))
    #  write(data, outfile, append=TRUE, ncolumn=length(data), sep="\t")
    #    }
    #    write_order_of_dysfunction('Order_of_dysfunction.txt', env, cells, isFirst)
    
    x <- array("",dim = c ( length(order_dysfunction[,1]) , 2) )
    
    for (i in 1:length(order_dysfunction[,1])) {
        data <- c(names(sort(order_dysfunction[i,2:(length(onco$name)+1)])))
        x[i,1] <- paste(data,collapse = " -> ")
        x[i,2] <- order_dysfunction$Num_cells[i]
    }
    
    x <- as.data.frame(x)
    x[,2] <- as.numeric(as.character(x[,2]))
    
    
    # print("The order of gene dysfunction for each clone is saved into the file")
    print("This is function to save the order of gene dysfunction to the file `Order_of_dysfunction.txt` ")
    # The order of gene dysfunction 
    
    # find the unique orders of genes dysfunction
    uniq_order <<- table(x[,1])
    uniq_order <<- sort(uniq_order,decreasing = TRUE)
    
    cells_order <<- vector(mode = "integer", length =  length(uniq_order))
    
    print("Unique of order of genes dysfunction and it's frequency (the first is a number of clones,  the second is number of cells):")
    for (i in 1:length(uniq_order)) {
        cells_order[i] <<- sum( x[ which(x[,1] == names(uniq_order)[i] ), 2] )
        
        print(c(names(uniq_order)[i], uniq_order[[i]], cells_order[i]))
        data <- c(names(uniq_order)[i], uniq_order[[i]], cells_order[i])
        write(data, outfile, append=TRUE, ncolumn=length(data), sep="\t")
    }
}  # End of function




######### Function to calulate TREE of CLUSTERS:

#### To plot tree of clusters
calc_tree <- function(evolution_clusters, clusters, total_clusters) { 
    cluster_tree <<- data.frame(matrix(nrow = total_clusters, ncol = 5))
    cluster_tree[,] <<- 0
    names(cluster_tree) <<- c("Cluster_ID","Time_start","Time_end","Parent_ID","Length")
    cluster_tree[,1] <<- names(evolution_clusters)      # cluster_ID
      
    for ( i in 1:total_clusters ) {
        st <- min(which(evolution_clusters[,i] > 0) )
        en <- max(which(evolution_clusters[,i] > 0) )
        cluster_tree[i,2] <<- st
        cluster_tree[i,3] <<- en    
        

        # To find parent ID: find first cell of cluster -> find parent ID of it -> 
        ##################################        ->   find cluster ID of parent cell at timestep  = birthday
        ### find first cell of cluster:
        ID <- as.numeric(cluster_tree[i,1])
        jk <- min(which(clusters[,"Cluster_ID"] == ID))
        P_ID <- clusters[jk,"Parent_ID"]   # cell ID of parent
        BirthDay <- clusters[jk,"Birthday"]  # birthday time step
        
        ## find cluster ID of parent cell at timestep  = birthday
        if (P_ID == 0)  { k <- min(which(clusters[,"Time"] == BirthDay) ) 
               } else {
                  k <- min(which( (clusters[,"Time"] == BirthDay) & 
                                    (clusters[,"Parent_ID"] == P_ID)))  # cell ID of parent at timestep
          }
        cluster_tree[i,"Parent_ID"] <<- clusters[k,"Cluster_ID"] 
    }
      
    cluster_tree[is.na(cluster_tree)] <<- 0 
      
    cluster_tree$Cluster_ID <<- as.numeric(cluster_tree$Cluster_ID)
    cluster_tree$Length <<- cluster_tree$Time_end - cluster_tree$Time_start
      
    cluster_tree$Length <<- cluster_tree$Length / max(cluster_tree$Length)
    
    cluster_tree$Length[ which(cluster_tree$Length == 0) ] <<- 0.01
      
    ### order(cluster_tree$Time_start) gives the order of appearence of clusters:
    cluster_tree[,] <<- cluster_tree[order(cluster_tree$Time_start),]
    cluster_tree[which(cluster_tree$Cluster_ID == 0),"Time_start"] <<- 0
      
    library(ape)
    n <- length(table(cluster_tree$Parent_ID))
      
    parents <- unlist(dimnames(table(cluster_tree$Parent_ID)))
    parents <- as.numeric(parents)
      
    children <- as.numeric(cluster_tree$Cluster_ID)
      
    m <- length(children)
      
    cl_edge <- data.frame(p=1:m,ch=1:m)
    cl_edge[,] <- -1
      
    k <- 1
    l_edge <- NULL
      
    for (i in 1:n) {
        vec_p <-  which(cluster_tree$Parent_ID == parents[i])
        for (j in vec_p) {
              
            if ( children[j] %in% parents ) {
                if ( children[j] != parents[1] )  {
                    cl_edge[k,] <- c(m+i, ( m + which(children[j] == parents)) )  
                    l_edge[k] <- 1.0
                    k <- k + 1
                }
                cl_edge[k,] <- c( (m + which(children[j] == parents)) , j )
                l_edge[k] <- 0.5
            }  else  { 
                cl_edge[k,] <- c(m+i,j) 
                l_edge[k] <- 0.5
            }
            k <- k + 1
        } 
    }
      
    tree_cl <- rtree(n = 2, rooted = TRUE)
      
    tree_cl$Nnode <- n 
    tree_cl$tip.label <- as.character( cluster_tree$Cluster_ID )
    tree_cl$edge <- as.matrix(cl_edge)
    tree_cl$edge.length <- l_edge
    attr(tree_cl$edge,"dimnames") <- NULL
      
      

    write.tree(tree_cl,file = "null")
      
    tree_cl <- read.tree(file = "null")
    tree_cl$root.edge <- 0.15
    file.remove("null")
      
    return(tree_cl)
}
    

################### Function to calculate the tree of CLONES:


#### To plot tree of clones
calc_tree_clones <- function(clones, total_clones) { 
  
    clones_new <- clones[which(clones[ ,"Clone_ID" ] < (total_clones +1)   ),]
    clone_tree <<- data.frame(matrix(nrow = total_clones, ncol = 5))
    clone_tree[,] <<- 0
    names(clone_tree) <<- c("Clone_ID","Time_start","Time_end","Parent_ID","Length")
    clone_tree[,1] <<- names(table(clones_new[,"Clone_ID"]))      # Clone_ID
    
    
    for ( i in 1:total_clones ) {
        st <- min( which(clones_new[,"Clone_ID"] == clone_tree[,"Clone_ID"][i]))
        en <- max( which(clones_new[,"Clone_ID"] == clone_tree[,"Clone_ID"][i]))
        clone_tree[i,"Time_start"] <<- clones_new[st,"Birthday"]
        clone_tree[i,"Time_end"]   <<- clones_new[en,"Time"]   
        clone_tree[i,"Parent_ID"]  <<- clones_new[st,"Parent_ID"] 
    }
    
    clone_tree[is.na(clone_tree)] <<- 0 
    
    clone_tree$Clone_ID <<- as.integer(clone_tree$Clone_ID)
    clone_tree$Length <<- clone_tree$Time_end - clone_tree$Time_start
    
    clone_tree$Length <<- clone_tree$Length / max(clone_tree$Length)
    
    clone_tree$Length[ which(clone_tree$Length == 0) ] <<- 0.01
    
    ### order(clone_tree$Time_start) gives the order of appearence of clones:
    clone_tree[,] <<- clone_tree[order(clone_tree$Time_start),]
    
    library(ape)
    n <- length(table(clone_tree$Parent_ID))
    
    parents <- unlist(dimnames(table(clone_tree$Parent_ID)))
    parents <- as.numeric(parents)
    
    children <- as.numeric(clone_tree$Clone_ID)
    
    m <- length(children)
    
    cl_edge <- data.frame(p=1:m,ch=1:m)
    cl_edge[,] <- -1
    
    k <- 1
    l_edge <- NULL
    
    for (i in 1:n) {
        vec_p <-  which(clone_tree$Parent_ID == parents[i])
        for (j in vec_p) {
            if ( children[j] %in% parents ) {
                if ( children[j] != parents[1] )  {
                    cl_edge[k,] <- c(m+i, ( m + which(children[j] == parents)) )  
                    l_edge[k] <- 1.0
                    k <- k + 1
                }
                cl_edge[k,] <- c( (m + which(children[j] == parents)) , j )
                l_edge[k] <- 0.5
              }  else  { 
                      cl_edge[k,] <- c(m+i,j) 
                      l_edge[k] <- 0.5
            }
            k <- k + 1
        } 
    }
    
    tree_cl <- rtree(n = 2, rooted = TRUE)
    
    tree_cl$Nnode <- n 
    tree_cl$tip.label <- as.character( clone_tree$Clone_ID )
    tree_cl$edge <- as.matrix(cl_edge)
    tree_cl$edge.length <- l_edge
    attr(tree_cl$edge,"dimnames") <- NULL
    
    
    
    write.tree(tree_cl,file = "null")
    
    tree_cl <- read.tree(file = "null")
    tree_cl$root.edge <- 0.15
    file.remove("null")
    
    return(tree_cl)
}




########### Inequality coefficient ############## Gini ################

# Functions to calculate the inequality coefficients:

ineq <- function(x, parameter = NULL, type=c("Gini", "RS", "Atkinson", "Theil",
                                             "Kolm", "var", "square.var", "entropy"), na.rm = TRUE)
{
  switch(match.arg(type),
         Gini = Gini(x, na.rm = na.rm),
         RS = RS(x, na.rm = na.rm),
         Atkinson = Atkinson(x, parameter = parameter, na.rm = na.rm),
         Theil = Theil(x, parameter = parameter, na.rm = na.rm),
         Kolm = Kolm(x, parameter = parameter, na.rm = na.rm),
         var = var.coeff(x, na.rm = na.rm),
         square.var = var.coeff(x, square=TRUE, na.rm = na.rm),
         entropy = entropy(x, parameter = parameter, na.rm = na.rm))
}

Gini <- function(x, corr = FALSE, na.rm = TRUE)
{
    if(!na.rm && any(is.na(x))) return(NA_real_)
    x <- as.numeric(na.omit(x))
    n <- length(x)
    x <- sort(x)
    G <- sum(x * 1L:n)
    G <- 2 * G/sum(x) - (n + 1L)
    if (corr) G/(n - 1L) else G/n
}

RS <- function(x, na.rm = TRUE)
{
    if(!na.rm && any(is.na(x))) return(NA_real_)
    x <- as.numeric(na.omit(x))
    d <- abs(x - mean(x))
    d <- mean(d)/(2*mean(x))
    d
}

Atkinson <- function(x, parameter = 0.5, na.rm = TRUE)
{
    if(!na.rm && any(is.na(x))) return(NA_real_)
    x <- as.numeric(na.omit(x))
    if(is.null(parameter)) parameter <- 0.5
    if(parameter==1)
      A <- 1 - (exp(mean(log(x)))/mean(x))
    else
    {
        x <- (x/mean(x))^(1-parameter)
        A <- 1 - mean(x)^(1/(1-parameter))
    }
    A
}

var.coeff <- function(x, square=FALSE, na.rm = TRUE)
{
  if(!na.rm && any(is.na(x))) return(NA_real_)
  x <- as.numeric(na.omit(x))
  n <- length(x)
  V <- sqrt((n-1)*var(x)/n)/mean(x)
  if(square) V <- V^2
  V
}

Theil <- function(x, parameter = 0, na.rm = TRUE)
{
    if(!na.rm && any(is.na(x))) return(NA_real_)
    x <- as.numeric(na.omit(x))
    if(is.null(parameter)) parameter <- 0
    if(parameter==0)
    {
        x <- x[!(x==0)]
        Th <- x/mean(x)
        Th <- sum(x*log(Th))
        Th <- Th/sum(x)
    }
    else
    {
        Th <- exp(mean(log(x)))/mean(x)
        Th <- -log(Th)
    }
    Th
}

Kolm <- function(x, parameter = 1, na.rm = TRUE)
{
    if(!na.rm && any(is.na(x))) return(NA_real_)
    x <- as.numeric(na.omit(x))
    x <- as.numeric(x)
    if(is.null(parameter)) parameter <- 1
    KM <- parameter * (mean(x)-x)
    KM <- mean(exp(KM))
    KM <- (1/parameter)*log(KM)
    KM
}

entropy <- function(x, parameter = 0.5, na.rm = TRUE)
{
    if(!na.rm && any(is.na(x))) return(NA_real_)
    x <- as.numeric(na.omit(x))
    x <- as.numeric(x)
    if(is.null(parameter)) parameter <- 0.5
    if(parameter==0)
      e <- Theil(x, parameter = 1)
    else
      if(parameter==1)
        e <- Theil(x, parameter = 0)
    else
    {
      k <- parameter
      e <- (x/mean(x))^k
      e <- mean(e-1)/(k*(k-1))
    }
    e
}


