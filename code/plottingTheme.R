# plotting theme
#######

cor_colors <- function(){

  colorVec <- c("#ff8080", "#1250cc", "#8c8c8c")
  names(colorVec) <- c("negRed","posBlue","zeroGray")
  return(colorVec)
}

blueGradient_colors <- function(){
  
  colorVec <- c("#1250cc", "#b7caef")
  names(colorVec) <- c("posBlue","low")
  return(colorVec)
}

rcp_colors <- function(){
  
  #"brown","orange","green","blue","pink"
  colorVec <- c("#6f2708","#dd7914","#1e632d","#072f5c","#b832c1")
  #old pink = #b256b8
  return(colorVec)
}

trophic.mode_colors <- function(){
  
  # The palette with black:
  # black, orange, lightblue, green, darkblue
  colorVec <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#0072B2")
  names(colorVec) <- c("unclassified", "Saprotroph", "Pathotroph","Symbiotroph","Pathotroph-Symbiotroph")
  
  return(colorVec)
}

make_ggplot_theme <- function(){
  
  require(ggplot2)
  
  #my ggplot template
  mytheme <- theme_bw(base_size = 10, base_family = "Helvetica") +
    theme(panel.border = element_rect(colour = "black"),      #put a black box around the plotting area
          axis.line = element_line(colour = "black"),                 #axis lines are in black
          panel.grid.major = element_blank(),                         #turn off the gridlines
          panel.grid.minor = element_blank(),
          strip.text.x = element_text(face='bold.italic', hjust=0.05),         #turn off the x axis facet labels
          strip.text.y = element_text(face='bold.italic', hjust=0.05),
          strip.background = element_rect(fill = 'white', colour='black'),    #make y axis facet labels be italic and top justified
          legend.key = element_blank(),                               #turn off box around legend
          plot.title=element_text(hjust=0, vjust=0.5, face='bold'), #style and position of the panel label
          plot.margin = unit(c(0.05,0.05,0.05,0.05),"in")
    )
  
  return(mytheme)
}

g_legend<-function(a.gplot){ 
  
  require(ggplot2)
  
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  
  return(legend)
} 


