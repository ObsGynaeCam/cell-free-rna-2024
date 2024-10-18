library(grDevices)
library(RColorBrewer)
source("libs/theme_publish.R") # defines theme_Publication()
source("libs/multiplot.R") # 'multiplot' is called within 'plotExp' below
							# called by 'bin/R/MJ.RNA/get.expression.profile.R'
							# ddsFpm should have been defined
############
## colour ##
############
hmcol.RdBu <- colorRampPalette(rev(brewer.pal(9, "RdBu")))(255)
hmcol.GnBu <- colorRampPalette(rev(brewer.pal(9, "GnBu")))(255)
hmcol<-hmcol.RdBu

# A colorblind-friendly palette
# http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
# 12 colours
# The palette with grey:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#662506", "#997a00", "#800080", "#000000")
#pie(rep(1, length(cbPalette)), col = cbPalette)
# The palette with black:
cbPalette2 <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
#pie(rep(1, length(cbPalette2)), col = cbPalette2)
#https://stackoverflow.com/questions/9563711/r-color-palettes-for-many-data-classes/41230685
c25 <- c(
            "black", "gold1",
            "skyblue2", "#FB9A99", # lt pink
            "palegreen2",
            "#FDBF6F", # lt orange
            "gray70", "khaki2",
            "#6A3D9A", # purple
            "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
            "darkturquoise", "green1", "yellow4", "yellow3",
            "darkorange4", "brown",
           "dodgerblue2", "#E31A1C", # red
            "green4",
            "#CAB2D6", # lt purple
            "#FF7F00" # orange
)
#pie(rep(1, 25), col = c25)

li.color=list(
    `source`=c('GTEx'=cbPalette[1],'Early gestation study'=cbPalette[2],'POP study'=cbPalette[3]),
    `GO`=c(`GO:BP`=ggsci::pal_d3()(3)[1],`GO:CC`=ggsci::pal_d3()(3)[2],`GO:MF`=ggsci::pal_d3()(3)[3]),
    `Any>0`=c(Yes="grey30",No="grey90"),
    #`Sex`=c(Female="#FB9A99",Male='skyblue2')
    `Sex`=c(Female="hotpink",Male='skyblue'),
    `Predicted Sex`=c(Female="hotpink3",Male='skyblue3')
)

############
# graphics #
############
#par(cex=0.5, cex.main=0.7, cex.lab=0.7)
#par(oma=c(0.2,0.2,0.2,0.2)) # outter margin area
#par(mar=c(5,5,2,2)) # margin area. the default is c(5, 4, 4, 2) + 0.1
myDotCol<-"#00000020" # transparent grey 

blankPlot <- ggplot()+geom_blank(aes(1,1))+
	theme(
		plot.background = element_blank(), 
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(), 
		panel.border = element_blank(),
		panel.background = element_blank(),
		axis.title.x = element_blank(),
		axis.title.y = element_blank(),
		axis.text.x = element_blank(), 
		axis.text.y = element_blank(),
		axis.ticks = element_blank(),
		axis.line = element_blank()
		)
