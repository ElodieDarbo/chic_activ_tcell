library(GenomicRanges)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(ggrepel)
library(igraph)
library(plotly)
#library(enrichR)
library(grr)
library(pheatmap)
library(Cairo) 

make.network.interaction <- function(interactions,RNAseq){
  temp <- interactions[,list(name1,name2,type1,type2,extend.classif,type)]
  net <- temp[,list(name1,name2,mode=extend.classif,type)]
  annots <- unique(as.data.frame(rbind(temp[,list(node=name1,type=type1)],temp[,list(node=name2,type=type2)])))
  G <- graph_from_data_frame(net)
  L <- as.data.frame(layout.kamada.kawai(G))
  
  m <- match(names(V(G)),annots$node)
  genes <- annots$node[na.omit(m)]

  V(G)$shape <- annots$type[match(genes,annots$node)]
  V(G)$label = annots$node[match(genes,annots$node)]
  V(G)$mode <- "node"
  vs <- V(G)
  es <- as.data.frame(get.edgelist(G))
  Nv <- length(vs)
  Ne <- nrow(es)
  
  L$name <- vs$label
  L$shape <- vs$shape
  
  m <- match(L$name,RNAseq$ID)
  L$mode <- "node"
  L$mode[!is.na(m)] <- RNAseq[na.omit(m)]$status
  edge_shapes <- rbindlist(lapply(1:Ne,function(i,es,edges.g,L){
    v0 <- as.vector(es[i,]$V1)
    v1 <- as.vector(es[i,]$V2)
    mode <- as.vector(edges.g$mode[i])
    type <- as.vector(edges.g$type[i])
    edge_shape <- data.table(V1=L$V1[L$name==v0],V2=L$V2[L$name==v0],xend=L$V1[L$name==v1],yend=L$V2[L$name==v1],
                             mode,type,name=mode)
    edge_shape
  },es=es,edges.g=E(G),L))
  print(edge_shapes)
  print(head(L))
  g <- ggplot(L,aes(x=V1,y=V2,label=name,color=mode)) + 
    geom_segment(data=edge_shapes,aes(xend=xend,yend=yend,linetype=type,color=mode)) +
    geom_point(aes(shape=shape,color=mode),size = 6,fill="white") +
    scale_color_manual(values=c(binding=colours()[54], stable="grey", up=colours()[54],down=colours()[124],gain="red",loss=colours()[124],catalysis=colours()[642],node="black",ptmod=colours()[466],inhibition=colours()[614],activation=colours()[144]),guide=FALSE) + 
    scale_linetype_manual(values=c(PPI="solid",bo="dashed",bb="solid")) +
    scale_shape_manual(values=c(b=23,o=21)) + theme_void() 
  g
}


multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


