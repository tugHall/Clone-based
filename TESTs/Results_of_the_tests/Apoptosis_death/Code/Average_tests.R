


  data_out <- read.csv("Output/cloneout.txt", sep="\t")
  data_out[is.na(data_out)] <- ""
  #View(data_out)      
 
  save_fig <- function(file_input) {
    tck_w <- 0.03
    file_input <- paste0(substr(file_input,1,nchar(file_input)-4),".jpg")
    #dev.copy2eps(file = file_input, height = 10, width = 10) 
    dev.copy(jpeg,file_input, width = 10, height = 10, units = 'in', res = 300)
    dev.off()
    # par(xpd=TRUE, cex.lab=2, lwd = 2, mar = c(5, 5, 5, 5), tcl = 0.5, cex.axis = 1.75,  mgp = c(3, 0.6, 0))
    par(cex.axis=1.3, cex.lab=1.6, cex=1.4, mar=c(5.1, 4.1, 2.1, 2.1), mgp=c(2.3, 0.3, 0),xpd=TRUE, tck = tck_w,font.lab = 2)
    
  }
  
  
  
 ave_plot <- function(data_out, output = "N", lg = "Number of cells", xl, yl, nm1 = "Fig_1.jpg", nm2 = "Fig_2.jpg") { 
      mT  <- max(data_out$Time)
      mRP <- max(data_out$Repeat_number)
       
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
      
      legend("top", c("100 simulations", "average value"), bty = "n", col = c("black","red"), 
                     lty = c(NA,1),pch = c(1,NA),pt.cex = 0.15,lwd = c(NA,2.4), horiz = TRUE)
      
      save_fig(paste0("Figures/",nm1) )
      rd <- readline(prompt= paste0("This is the simple plot for ", lg, " - Press Enter  "))
      
      library(MASS)
      
      df <- data.frame(x=data_out$Time,y=d_out)
      
      # colors
      library(RColorBrewer)
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
            for (kk in 1:(mT-1)) {
            
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

