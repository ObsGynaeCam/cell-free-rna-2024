# code from https://rpubs.com/Koundy/71792

#theme_Publication <- function(base_size=14, base_family="helvetica") {
theme_Publication <- function(base_size=14, base_family="") {
      library(grid)
      library(ggthemes)
      (theme_foundation(base_size=base_size, base_family=base_family)
       + theme(
            plot.title = element_text(face = "bold",size = rel(1.3), hjust = 0.5),
			text = element_text(),
			panel.background = element_rect(colour = NA),
			plot.background = element_rect(colour = NA),
			panel.border = element_rect(colour = NA),
			axis.title = element_text(face = "bold",size = rel(1.2)),
			axis.title.y = element_text(angle=90,vjust =2),
			axis.title.x = element_text(vjust = -0.2),
			axis.text = element_text(size=rel(1.1)), 
			axis.line = element_line(colour="black"),
			axis.ticks = element_line(),
			panel.grid.major = element_line(colour="#f0f0f0"),
			panel.grid.minor = element_blank(),
			legend.key = element_rect(colour = NA),
			#legend.background = element_rect(colour = 'black', fill='grey',linetype='dashed'),
			legend.position = "right",
			#legend.direction = "horizontal",
			legend.key.size= unit(0.7, "cm"),
			legend.spacing= unit(0.1, "mm"),
			legend.title = element_text(face="bold.italic",size=rel(1)),
			legend.text = element_text(size = rel(0.9),family = "sans"),
			plot.margin=unit(c(10,5,5,5),"mm"),
			#strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"), # grey
      strip.background = element_rect(fill = "#17252D", color = "#17252D"), # black-ish
			strip.text = element_text(face="bold",color="white", size=rel(1.1))
          ) + if(packageVersion("ggplot2")<=2.1){theme(legend.margin = unit(0.1, "mm"))}else{theme(legend.spacing = unit(0.1, "mm"))}
	   )
      
}

# https://benjaminlouis-stat.fr/en/blog/2020-05-21-astuces-ggplot-rmarkdown/
theme_Sung <- function(base_size = 14) {
  theme_bw(base_size = base_size) %+replace%
    theme(
      # Title
      plot.title = element_text(size = rel(1), face = "bold", margin = margin(0,0,5,0), hjust = 0),
      # Zone
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      # axes
      axis.title = element_text(size = rel(0.85), face = "bold"),
      axis.text = element_text(size = rel(0.70), face = "bold"),
      axis.line = element_line(color = "black", arrow = arrow(length = unit(0.3, "lines"), type = "closed")),
      # legend
      legend.title = element_text(size = rel(0.85), face = "bold"),
      legend.text = element_text(size = rel(0.70), face = "bold"),
      legend.key = element_rect(fill = "transparent", colour = NA),
      legend.key.size = unit(1.5, "lines"),
      legend.background = element_rect(fill = "transparent", colour = NA),
      # facetting
      strip.background = element_rect(fill = "#17252D", color = "#17252D"),
      strip.text = element_text(size = rel(0.85), face = "bold", color = "white", margin = margin(5,0,5,0))
    )
}

scale_fill_Publication <- function(...){
      library(scales)
      discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)

}

scale_colour_Publication <- function(...){
      library(scales)
      discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)

}

