
### create GIF animation:

dr <- '../Figures/TIME_CELL_AND_CELL/'

setwd(dr)

# example 1: simple animated countdown from 10 to "GO!".
ex <- FALSE

if (ex) {
    png(file="example%02d.png", width=200, height=200)
    for (i in c(10:1, "G0!")){
      plot.new()
      text(.5, .5, i, cex = 6)
    }
    dev.off()
}

# convert the .png files to one .gif file using ImageMagick. 
# The system() function executes the command as if it was done
# in the terminal. the -delay flag sets the time between showing
# the frames, i.e. the speed of the animation.

# install.packages("magick")
library("magick")

if (ex) {
    ims <- c("example01.png","example02.png",  
             "example03.png","example04.png", "example05.png",
             "example06.png","example07.png", "example08.png",
             "example09.png","example10.png", "example11.png") 
}

b_nm <- 'env_metastasis_' 
ims <- NULL

for (p in 1:9 ) {
    ims <- c(ims, paste0(b_nm,p,'.jpg') )
}


mg <- image_read( ims )


an <- image_animate(mg, fps = 10)

# image_animate(mg, fps = 10)


## install.packages("gifski")


library("gifski")

image_write_gif(an, path =paste0( b_nm, 'animate', '.gif') , delay = 1)




# to not leave the directory with the single jpeg files
# I remove them.


# file.remove(list.files(pattern=".png"))
