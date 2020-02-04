
### Make plot 
tck_w <- 0.03
par(cex.axis=1.3, cex.lab=1.6, cex=1.4, mar=c(5.1, 4.1, 2.1, 2.1), mgp=c(2.3, 0.3, 0),xpd=TRUE, tck = tck_w,font.lab = 2)

library(MASS)
library(RColorBrewer)

library(ggplot2)
library(plyr)


save_fig <- function(file_input) {
    #tck_w <- 0.03
    # file_input <- paste0(substr(file_input,1,nchar(file_input)-4),".jpg")
    #dev.copy2eps(file = file_input, height = 10, width = 10) 
    dev.copy(jpeg,file_input, width = 10, height = 10, units = 'in', res = 300)
    dev.off()
    # par(xpd=TRUE, cex.lab=2, lwd = 2, mar = c(5, 5, 5, 5), tcl = 0.5, cex.axis = 1.75,  mgp = c(3, 0.6, 0))
    # par(cex.axis=1.3, cex.lab=1.6, cex=1.4, mar=c(5.1, 4.1, 2.1, 2.1), mgp=c(2.3, 0.3, 0),xpd=TRUE, tck = tck_w,font.lab = 2)
  
}


ave_plot <- function(data_out, output = "N", lg = "Number of cells", xl, yl, nm1 = "Fig_1.jpg", nm2 = "Fig_2.jpg") { 
    mT  <- max(data_out$Time)
    minT <- min(data_out$Time)
    mRP <- max(data_out$Repeat_number)
    
    if (mT > xl[2] ) {mT <- xl[2] - 1}
    
    tck_w <- 0.03
    par(cex.axis=1.3, cex.lab=1.6, cex=1.4, mar=c(5.1, 4.1, 2.1, 2.1), mgp=c(2.3, 0.3, 0),xpd=TRUE, tck = tck_w,font.lab = 2)
    
    d_out <- data_out[[output]]
    # plot(data_out$Time,d_out,type = "p",cex = 0.2)
    if (  sum(d_out) == 0) { stop("There are all zero data, please, check input data !") }
    
    
    data_ave <- matrix(data = 0, ncol = 2,nrow = mT + 1)
    for (k in 1:(mT+1)  ) {
        data_ave[k,1] = k-1
        n <- which(data_out$Time == (k-1) )
        data_ave[k,2] = sum(d_out[n]) / mRP
    }
    
    
    plot(data_out$Time,d_out,type = "p",cex = 0.15, xlab = "Time step", ylab = lg, xlim = xl, ylim = yl)
    axis(3,tck = tck_w,labels = FALSE)
    axis(4,tck = tck_w,labels = FALSE)
    
    
    lines(data_ave[1:(mT+1),1],data_ave[1:(mT+1),2],type = "l",lwd = 3, col = "red", add = TRUE)
    
    legend("top", c('100 simulations        ', 'average value'), bty = "n", col = c("black","red"), 
           lty = c(NA,1),pch = c(1,NA),pt.cex = 0.15,lwd = c(NA,2.4), horiz = TRUE)
    
    save_fig(paste0("Figures/",nm1) )
    rd <- readline(prompt= paste0("This is the simple plot for ", lg, " - Press Enter  "))
    

    
    df <- data.frame(x=data_out$Time,y=d_out)
    
    # colors

    rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
    r <- rf(40)
    
    # Adjust binning (interpolate - can be computationally intensive for large datasets)
    k <- kde2d(df$x, df$y, n=c(300,300))
    k$z <- k$z / max(k$z)
    
    breaks <- 1:(length(r)+1)
    breaks <- breaks / length(r)
    #image(k, col=r, breaks = breaks, xlab = "Generation number", ylab = "Number of cells",ylim = c(0,1500))
    #axis(3,tck = tck_w,labels = FALSE)
    #axis(4,tck = tck_w,labels = FALSE)
    
    #lines(data_ave[2:100,1],data_ave[2:100,2],type = "l",lwd = 2.4, col = "#FF3300")
    
    #levels <- c(min(breaks),0.05,0.1,0.2,0.4,0.6,0.8)
    #contour(k, levels = levels, lty = 1 , lwd=0.9, col= "#00EE99", drawlabels = FALSE, add = TRUE)
    
    
    for (i in 1:2) {
        image.default(k, col=r, breaks = breaks, xlab = "Time step", ylab = lg, ylim = yl, xlim = xl, legend=TRUE)
        
        # For each time step
        for (kk in (minT+1):(mT-1)) {
            numb <- which( xor(df$x==kk,df$x==(kk+1)) )
            h = ifelse( max(df$y[numb]) - min(df$y[numb]) > 0, max(df$y[numb]) - min(df$y[numb]), 1)
            k <- kde2d(df$x[numb], df$y[numb], h = c(1,h), n=c(2,100))
            k$z <- k$z / max(k$z)
            
            breaks <- 1:(length(r)+1)
            breaks <- breaks / length(r)
            image(k, col=r, breaks = breaks, add = TRUE)
        }
    }      
    axis(3,tck = tck_w,labels = FALSE)
    axis(4,tck = tck_w,labels = FALSE)
    
    lines(data_ave[1:(mT+1),1],data_ave[1:(mT+1),2],type = "l",lty = 1, lwd = 3, col = "#FF00FF")
    
    legend(x = xl[2] * 0.6, y = yl[2] * 0.97, c("average value",""), bty = "n", col = c("#FF00FF",NA),lty = c(1,NA),pch = c(NA,NA),pt.cex = 0.15,lwd = c(4,NA), horiz = TRUE)
    save_fig(paste0("Figures/",nm2) )
    
}





### Time of simulations:
TM <- data.frame(cell_based = c(257,494,261,493,270), 
                 clone_based = c(5.5, 8.5, 5.2, 8.9, 5.45),
                 row.names = c("Apoptosis_primary", "Apoptosis_metastasis",
                               "Environment_primary", "Environment_metastasis",
                               "Metastasis_transformation") )

### cell-based is calculated on 9 nodes with 24 cores each, so we should x by 216  (9x24=216):

TM$cell_based <- TM$cell_based * 216


### Acceleration plot:

barplot(TM$cell_based/TM$clone_based, col = "blue", cex.names = 1.0, space = c(0.6,0.4),
        names.arg = c("Ap_Pr", "Ap_Met", "Env_Pr", "Env_Met", "Met_Tr"), 
        ylab = 'Acceleration rate') 

nm <- 'Acceleration.jpg'
save_fig(paste0("Figures/",nm) )


### I) Comparison for tests about Apoptosis death for metastasis cells:

DT_CLONES <- read.csv2(file = 'RESULTS_CLONES/apoptosis_clones_metastasis.txt', sep = '\t')
DT_CLONES[is.na(DT_CLONES)] <- ''        

ave_plot(data_out = DT_CLONES, output = 'M',lg = 'Number of metastasis cells', xl = c(0,100), yl = c(0,2500), 
         nm1 = 'Apoptosis_clones_metastasis_1.jpg', nm2 = 'Apoptosis_clones_metastasis_2.jpg')


DT_CELLS <- read.csv2(file = 'RESULTS_CELLS/apoptosis_cells_metastasis.txt', sep = ',')
DT_CELLS[is.na(DT_CELLS)] <- ''
names(DT_CELLS)[which(names(DT_CELLS) == 'trial')] <- 'Repeat_number'

ave_plot(data_out = DT_CELLS, output = 'M',lg = 'Number of metastasis cells', xl = c(0,100), yl = c(0,2500), 
         nm1 = 'Apoptosis_cells_metastasis_1.jpg', nm2 = 'Apoptosis_cells_metastasis_2.jpg')



### II) Comparison for tests about Apoptosis death for primary cells:

DT_CLONES <- read.csv2(file = 'RESULTS_CLONES/apoptosis_clones_primary.txt', sep = '\t')
DT_CLONES[is.na(DT_CLONES)] <- ''        

ave_plot(data_out = DT_CLONES, output = 'N',lg = 'Number of primary cells', xl = c(0,80), yl = c(0,2200), 
         nm1 = 'Apoptosis_clones_primary_1.jpg', nm2 = 'Apoptosis_clones_primary_2.jpg')


DT_CELLS <- read.csv2(file = 'RESULTS_CELLS/apoptosis_cells_primary.txt', sep = ',')
DT_CELLS[is.na(DT_CELLS)] <- ''
names(DT_CELLS)[which(names(DT_CELLS) == 'trial')] <- 'Repeat_number'

ave_plot(data_out = DT_CELLS, output = 'N',lg = 'Number of primary cells', xl = c(0,80), yl = c(0,2200), 
         nm1 = 'Apoptosis_cells_primary_1.jpg', nm2 = 'Apoptosis_cells_primary_2.jpg')



### III) Comparison for tests about Environment death for metastasis cells:

DT_CLONES <- read.csv2(file = 'RESULTS_CLONES/environment_clones_metastasis.txt', sep = '\t')
DT_CLONES[is.na(DT_CLONES)] <- ''        

ave_plot(data_out = DT_CLONES, output = 'M',lg = 'Number of metastasis cells', xl = c(0,100), yl = c(0,2500), 
         nm1 = 'Environment_clones_metastasis_1.jpg', nm2 = 'Environment_clones_metastasis_2.jpg')


DT_CELLS <- read.csv2(file = 'RESULTS_CELLS/environ_cells_metastasis.txt', sep = ',')
DT_CELLS[is.na(DT_CELLS)] <- ''
names(DT_CELLS)[which(names(DT_CELLS) == 'trial')] <- 'Repeat_number'

ave_plot(data_out = DT_CELLS, output = 'M',lg = 'Number of metastasis cells', xl = c(0,100), yl = c(0,2500), 
         nm1 = 'Environment_cells_metastasis_1.jpg', nm2 = 'Environment_cells_metastasis_2.jpg')


### IV) Comparison for tests about Environment death for primary cells:

DT_CLONES <- read.csv2(file = 'RESULTS_CLONES/environment_clones_primary.txt', sep = '\t')
DT_CLONES[is.na(DT_CLONES)] <- ''        

ave_plot(data_out = DT_CLONES, output = 'N',lg = 'Number of primary cells', xl = c(0,70), yl = c(0,2200), 
         nm1 = 'Environ_clones_primary_1.jpg', nm2 = 'Environ_clones_primary_2.jpg')


DT_CELLS <- read.csv2(file = 'RESULTS_CELLS/environ_cells_primary.txt', sep = ',')
DT_CELLS[is.na(DT_CELLS)] <- ''
names(DT_CELLS)[which(names(DT_CELLS) == 'trial')] <- 'Repeat_number'

ave_plot(data_out = DT_CELLS, output = 'N',lg = 'Number of primary cells', xl = c(0,70), yl = c(0,2200), 
         nm1 = 'Environ_cells_primary_1.jpg', nm2 = 'Environ_cells_primary_2.jpg')


### V) Comparison for tests about Invasion/metastasis transformation:

DT_CLONES <- read.csv2(file = 'RESULTS_CLONES/meta_clones_transformation.txt', sep = '\t')
DT_CLONES[is.na(DT_CLONES)] <- ''        

ave_plot(data_out = DT_CLONES, output = 'N',lg = 'Number of primary cells', xl = c(0,70), yl = c(0,2200), 
         nm1 = 'Meta_transform_clones_1.jpg', nm2 = 'Meta_transform_clones_2.jpg')


DT_CELLS <- read.csv2(file = 'RESULTS_CELLS/meta_cells_transformation.txt', sep = ',')
DT_CELLS[is.na(DT_CELLS)] <- ''
names(DT_CELLS)[which(names(DT_CELLS) == 'trial')] <- 'Repeat_number'

ave_plot(data_out = DT_CELLS, output = 'N',lg = 'Number of primary cells', xl = c(0,70), yl = c(0,2200), 
         nm1 = 'Meta_transform_cells_1.jpg', nm2 = 'Meta_transform_cells_2.jpg')






##############################################################################################################
tck_w <- 0.03
par(cex.axis=1.3, cex.lab=1.6, cex=1.4, mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(2.3, 0.3, 0),xpd=TRUE, tck = tck_w,font.lab = 2)

my_plot <- function(DT, xlab = 'x', ylab = 'y') {
  ### DT is dataframe with data and source columns: DT$data , DT$source 
  cdat <- ddply(DT, "source", summarise, data.mean=mean(data))
  
       bck <- theme_bw() + 
              theme(
              axis.text = element_text(face="bold", color="#000000", size=14, angle=0),
              axis.title = element_text(face="bold", color="#000000", size=18, angle=0),
              axis.line = element_line(colour = "darkblue", margin(2,2,2,2,unit = 'cm'),
                                       size = 1, linetype = "solid"),
#             panel.background = element_rect(fill = 'white', colour = 'black' )
#              axis.ticks.x.top =element_line(color = 'black', size = 0.7 )
#              panel.border = element_rect(fill = NULL, colour = '#0000FF', size = 2 )
#              axis.ticks = element_line(color = 'black', size = 0.7 ),
              axis.ticks.length=unit(-0.25, "cm"),
              axis.text.x = element_text(margin = margin(t = .5, unit = "cm")),
              axis.text.y = element_text(margin = margin(r = .5, unit = "cm")),
              legend.position = c(.95, .95),
              legend.justification = c("right", "top"),
              legend.text = element_text(size = 12, colour = "black"),
              legend.title = element_text(size = 18, color = 'red', face = 'bold'),
              plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
#              axis.ticks.margin = unit(2.5, "cm") 
    ) 

    ggplot(DT, aes(x=data, colour=source)) +
        geom_density(size = 1.6) +
        geom_vline(data=cdat, aes(xintercept=data.mean,  colour=source),
                   linetype="dashed", size=0.8) + 
        xlab(xlab) + ylab(ylab) + bck + 
      scale_x_continuous(sec.axis = dup_axis(labels = NULL, name = NULL) ) +
      scale_y_continuous(sec.axis = dup_axis(labels = NULL, name = NULL) )
}


### Compare distributions: 
shift <- 2
output <- 'M'
dr <- paste0("Figures/",'TIME/')    # directory to save the series of figures in to 'Figures' folder
nm <- 'env_metastasis_'  # the base of name for series of figures

if ( !dir.exists(dr) ) { dir.create(dr)  }

for (p in 1:9) {

    t <- 10 * p
    CE <- DT_CELLS[which(DT_CELLS$Time == t + shift & DT_CELLS$Repeat_number < 541),]
    # CL <- DT_CLONES[which(DT_CLONES$Time == t & DT_CLONES$Repeat_number < 1081),]
    CL <- DT_CELLS[which(DT_CELLS$Time == t + shift & DT_CELLS$Repeat_number > 540),]
    
    DT1 <- data.frame(data = CE[[output]])
    DT1$source <- 'CELL-based < 540'
    DT2 <- data.frame(data = CL[[output]])
    # DT2$source <- 'CLONE-based'
    DT2$source <- 'CELL-based > 540'
    DT <- rbind(DT1,DT2)
    
    print( my_plot(DT, xlab = 'Number of cells', ylab = 'Distribution') )
    
    save_fig(paste0(dr,nm,p,'.jpg') )
    rd <- readline(prompt= paste0("This is the time series of distributions for cell-based and clone-based data. Time is  ", t, " - Press for next  "))
    
}





