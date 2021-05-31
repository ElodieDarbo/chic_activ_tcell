main.dir <- ""
#main.dir <- "/Users/elodiedarbo/Documents/projects/C-HiC_cell_T/app_integration"
res.dir <- data.dir <- file.path(main.dir,"data")

options(shiny.maxRequestSize=50*1024^2) 

shinyServer(function(input, output) {
  ########################################################################
  ## DATA
  ########################################################################
    withProgress(message = 'loading data', {
      message("Loading data ...")
      source(file.path(main.dir,"scripts","functions_for_report.R"))
      load(file.path(data.dir,"refseq_hg19.RData"))
      # homer.tmp,shaps,tsne.group,tsne.center,peaks.ATAC.frag,all.interactions,
      # course.expression,homer.num,external.chip.ATAC,YY1,TF,peaks.ATAC,RNAseq,tad_hg19, 
      # tads_borders
      load(file.path(res.dir,"objects_for_integration.RData")) 
      time.course <- data.table(read.table(file.path(data.dir,"time_course_GSE138767.tab"),sep="\t",head=T))
      files.H3 <- system(paste("ls",file.path(data.dir,"ChIP_seq","*summary")),intern=T)
      H3 <- do.call("c",lapply(files.H3,function(f){
        nom <- unlist(sapply(strsplit(basename(f),"-"),"[[",1))
        cond <- unlist(sapply(strsplit(nom,"_"),"[[",1))
        mark <- unlist(sapply(strsplit(nom,"_"),"[[",2))
        tmp <- read.table(f,sep="\t",head=T)
        colnames(tmp) <- c("seqnames","start","end","nb.cond","nb.control","pval","score","FDR")
        tmp$cond <- cond
        tmp$mark <- mark
        tmp <- as(tmp,"GRanges")
        tmp <- tmp[tmp$FDR<0.01]
      }))
      load(file.path(data.dir,"YY1_CNR.RData"))
      peaks.ATAC$H3K27me3 <- peaks.ATAC$H3K4me1 <- peaks.ATAC$H3K4me3 <- peaks.ATAC$H3K27Ac <- NULL
      lapply(c("H3K27me3","H3K4me1","H3K4me3","H3K27Ac"),function(mark,H3){
        G0 <- H3[H3$mark==mark & H3$cond=="G0"]
        G1 <- H3[H3$mark==mark & H3$cond=="G1"]
        ovl.G0 <- findOverlaps(peaks.ATAC,G0,minoverlap = 25)
        ovl.G1 <- findOverlaps(peaks.ATAC,G1,minoverlap = 25)
        pos <- 1:length(peaks.ATAC)
        posG1 <- queryHits(ovl.G1)
        posG0 <- queryHits(ovl.G0)
        peaks.ATAC$presence <<- ifelse(!pos%in%c(posG1,posG0),"00",ifelse(pos%in%posG1[posG1%in%posG0],"11",ifelse(pos%in%posG1[!posG1%in%posG0],"01","10")))
        colnames(mcols(peaks.ATAC))[ncol(mcols(peaks.ATAC))] <<- mark
        mark
      },H3)
      
      load(file.path(data.dir,"external.datasets.RData"))
      extern.chip <- c(GSE62486.peaks,GSE116695.peaks,list(YY1.down=YY1[YY1$diff=="down"],YY1.up=YY1[YY1$diff=="up"]))
      annots.extern.chip <- rbind(data.table(GSE62486.annots)[,list(ID,TF,rep=NA,cell,treatment)],data.table(GSE116695.annots)[,list(ID,TF,rep,cell=NA,treatment=NA)],data.table(ID=c("YY1.down","YY1.up"),TF=c("YY1","YY1"),rep=NA,cell=NA,treatment=NA))
      
      p50 <- reduce(c(GSE116695.peaks$p50_1,GSE116695.peaks$p50_3))
      NFAT1 <- reduce(c(GSE116695.peaks$NFAT1_1,GSE116695.peaks$NFAT1_3))
      NFAT2 <- reduce(c(GSE116695.peaks$NFAT2_1,GSE116695.peaks$NFAT2_3))
      polII_RA_Th1 <- reduce(c(GSE62486.peaks$Th1_reactivated_PolII,GSE62486.peaks$Th1_RA_PolII))
      polII_RA_Th2 <- reduce(c(GSE62486.peaks$Th2_reactivated_PolII,GSE62486.peaks$Th2_RA_PolII))
      PolII_US_Th1 <- reduce(c(GSE62486.peaks$Th1_unstimulated_PolII,GSE62486.peaks$Th1_US_PolII))
      PolII_US_Th2 <- reduce(c(GSE62486.peaks$Th2_unstimulated_PolII,GSE62486.peaks$Th2_US_PolII))
      
      exclude.ds <- c("Th1_Native_H3K4me3","Th1_RA_IgG","Th2_RA_IgG","Th2_Native_H3K4me3","p50_1","p50_3","NFAT1_1","NFAT2_3","NFAT2_1","NFAT1_3","Th1_reactivated_PolII","Th1_RA_PolII","Th2_reactivated_PolII","Th2_RA_PolII","Th1_unstimulated_PolII","Th1_US_PolII","Th2_unstimulated_PolII","Th2_US_PolII")
      
      extern.chip <- extern.chip[!names(extern.chip)%in%exclude.ds]
      annots.extern.chip <- annots.extern.chip[!ID%in%exclude.ds]
      
      extern.chip <- c(extern.chip,p50,NFAT1,NFAT2,polII_RA_Th1,polII_RA_Th2,PolII_US_Th1,PolII_US_Th2)
      annots.extern.chip <- rbind(annots.extern.chip,data.table(ID= c("p50","NFAT1","NFAT2","polII_RA_Th1","polII_RA_Th2","PolII_US_Th1","PolII_US_Th2"),TF=NA,rep=NA, cell=NA,treatment=NA))
      
      external.chip.ATAC <- lapply(extern.chip,function(x,peaks.ATAC){
        ovl <- findOverlaps(peaks.ATAC,x,minoverlap = 25,select="first")
        res <- data.table(chip=ifelse(is.na(ovl),0,1))
        colnames(res) <- ""
        res
      },peaks.ATAC)
      external.chip.ATAC <- do.call("cbind",external.chip.ATAC)
      colnames(external.chip.ATAC) <- annots.extern.chip$ID
      external.chip.ATAC$ATAC.id <- peaks.ATAC$id
      message("Loading data done")
      
    })
  
  #################################################################################################
  ## accueil
  #################################################################################################
  
  output$acceuilText <- renderText({
    txt <- "Application is ready to use."
  })
  
  
  
  #################################################################################################
  ## whole data 
  #################################################################################################
  
  
  tab.interaction <- reactive({
    message("Preparing table with interaction ...")
    interactions <- all.interactions[,list(name1,id1,name2,id2,diff=extend.classif,B2B=type,bait1=RNAseq[match(name1,ID)]$status,bait2=RNAseq[match(name2,ID)]$status)]
    interactions$bait1[is.na(interactions$bait1)] <- "not.in.data"
    interactions$bait2[is.na(interactions$bait2) & all.interactions$type2=="b"] <- "not.in.data"
    interactions$distance <- all.interactions$distance
    if (input$only.diff){
      interactions <- interactions[diff!="stable"]
    }
    interactions
  })
  
  output$tab.interaction <- renderDataTable({
    tab.interaction()
  })
  
  write.tab.interaction <- eventReactive(input$save.tab.int,{
    tab <- tab.interaction()
    if (input$only.diff){
      if (file.exists(file.path(data.dir,"tables","differential.interactions.tab"))){
        text <- "differential.interactions.tab already exists in the tables folder"
      }
      else {
        write.table(tab,file.path(data.dir,"tables","differential.interactions.tab"),sep="\t",quote=F,row.names=F)
        text <- "Interactions have been saved in differential.interactions.tab"
      }
    }
    else {
      if (file.exists(file.path(data.dir,"tables","all.interactions.tab"))){
        text <- "all.interactions.tab already exists in the tables folder"
      }
      else {
        write.table(tab,file.path(data.dir,"tables","all.interactions.tab"),sep="\t",quote=F,row.names=F)
        text <- "Interactions have been saved in all.interactions.tab"
      }
    }
    text
  })
  
  output$write.tab.interaction <- renderText({
    write.tab.interaction()
  })
  
  per.fragment.interactions.total <- reactive({
    message("Preparing table with fragments ...")
    int.fragments.stats <- int.fragments()
    int.fragments.stats <- int.fragments.stats[,list(total=length(diff.int),
                                                     nb.gained=sum(diff.int=="gain"),
                                                     nb.lost=sum(diff.int=="loss"),
                                                     nb.differential=sum(diff.int!="stable"),
                                                     direction=sum(direction),
                                                     direction.diff=sum(direction[diff.int!="stable"]),
                                                     directionality = round(sum(direction)/length(diff.int),2)
    ),by=c("name","type")]
    int.fragments.stats$chr <- "0"
    int.fragments.stats$start <- 0
    int.fragments.stats$pos <- "0"
    int.fragments.stats[type=="b"]$chr <- as.vector(seqnames(all.fragments))[match(int.fragments.stats[type=="b"]$name,all.fragments$gene)]
    int.fragments.stats[type=="b"]$start <- as.vector(start(all.fragments))[match(int.fragments.stats[type=="b"]$name,all.fragments$gene)]
    int.fragments.stats[type=="b"]$pos <- as.vector(all.fragments$pos)[match(int.fragments.stats[type=="b"]$name,all.fragments$gene)]
    int.fragments.stats
  })
  
  tab.fragments <- reactive({
    tab <- per.fragment.interactions.total()
    tab$pos <- NULL
    if (input$only.diff.frag){
      tab <- tab[nb.differential>0]
    }
    tab
  })
  
  output$tab.fragments <- renderDataTable({
    tab.fragments()
  })
  
  write.tab.fragments <- eventReactive(input$save.tab.frag,{
    tab <- tab.fragments()
    if (input$only.diff.frag){
      if (file.exists(file.path(data.dir,"tables","frag_with_diff_interaction.tab"))){
        text <- "frag_with_diff_interaction.tab already exists in the tables folder"
      }
      else {
        write.table(tab,file.path(data.dir,"tables","frag_with_diff_interaction.tab"),sep="\t",quote=F,row.names=F)
        text <- "Interactions have been saved in frag_with_diff_interaction.tab"
      }
    }
    else {
      if (file.exists(file.path(data.dir,"tables","frag_with_interaction.tab"))){
        text <- "frag_with_interaction.tab already exists in the tables folder"
      }
      else {
        write.table(tab,file.path(data.dir,"tables","frag_with_interaction.tab"),sep="\t",quote=F,row.names=F)
        text <- "Interactions have been saved in frag_with_interaction.tab"
      }
    }
    text
  })
  
  output$write.tab.fragments <- renderText({
    write.tab.fragments()
  })
  
  
  all.ATAC.peaks.annots <- reactive({
    message("Preparing table with ATAC peaks ...")
    ATAC <- data.table(as.data.frame(peaks.ATAC))[,list(ID=id,status,genomic,SYMBOL)]
    ovl <- findOverlaps(peaks.ATAC,all.fragments)
    ovl <- data.table(as.data.frame(ovl))
    frags <- all.fragments[ovl$subjectHits]
    ovl$frag.id <- ifelse(frags$is.bait=="other",frags$IDs,frags$gene)
    ovl <- ovl[,list(frag.id=paste(unique(frag.id),collapse=", ")),by=queryHits]
    ATAC$frag.id[ovl$queryHits] <- ovl$frag.id
    #homer.tmp <-  prepare.shaps()$homer.tmp
    homer.tmp$db.clust <- as.vector(homer.tmp$db.clust)
    homer.tmp$db.clust[!is.na(homer.tmp$cl1.split)] <- homer.tmp$cl1.split[!is.na(homer.tmp$cl1.split)] 
    m <- data.table(as.data.frame(matches(ATAC$ID,homer.tmp$ATAC.id)))
    m$db.clust <- homer.tmp$db.clust[m$y]
    m <- m[,list(db.clust=paste(unique(db.clust),collapse=", ")),by=x]
    m$db.clust[is.na(m$db.clust) | m$db.clust=="NA"] <- ""
    ATAC$tSNE[m$x] <- m$db.clust
    ATAC <- cbind(ATAC,as.data.frame(peaks.ATAC)[,grep("H3",colnames(as.data.frame(peaks.ATAC)))])
    external.chip.ATAC <- as.data.frame(external.chip.ATAC)
    chip <- colnames(external.chip.ATAC)[-ncol(external.chip.ATAC)]
    ATAC$CHiP <- apply(external.chip.ATAC,1,function(x,chip){
      paste(unique(chip[x>0]),collapse=", ")
    },chip)
    colnames(homer.num) <- unlist(sapply(strsplit(colnames(homer.num),"[.][.]"),"[[",1))
    homer.num <- homer.num[match(ATAC$ID,homer.num$ATAC.id),]
    TFs <- colnames(homer.num)[1:436]
    ATAC$TFBSs <- apply(homer.num,1,function(x,TFs){
      paste(unique(TFs[x>0]),collapse=", ")
    },TFs)
    ATAC
  })
  
  output$tab.ATAC <- renderDataTable({
    all.ATAC.peaks.annots()
  })
  
  write.tab.ATAC <- eventReactive(input$save.tab.ATAC,{
    tab <- all.ATAC.peaks.annots()
    if (file.exists(file.path(data.dir,"tables","ATAC_detailed.tab"))){
        text <- "ATAC_detailed.tab already exists in the tables folder"
    }
    else {
        write.table(tab,file.path(data.dir,"tables","ATAC_detailed.tab"),sep="\t",quote=F,row.names=F)
        text <- "Interactions have been saved in ATAC_detailed.tab"
    }
    
    text
  })
  
  output$write.tab.ATAC <- renderText({
    write.tab.ATAC()
  })
  
  tab.expression <- reactive({
    colnames(time.course)[c(1,2)] <- c("ensEMBLID","ID")
    time.course <- time.course[!is.na(ID) & ID!=""]
    expression <- merge(RNAseq[,list(ID,status.RNAseq=status)],time.course[,list(ID,ensEMBLID,time.course,T0.T20,T20.T1,T1.T2,T2.T4,T4.T24)],by="ID",all=T)
    expression <- expression[,list(ID,status.RNAseq,time.course,T0.T20,T20.T1,T1.T2,T2.T4,T4.T24)]
  })
  
  output$tab.expression <- renderDataTable({
    tab.expression()
  })
  
  write.tab.expression <- eventReactive(input$save.tab.expr,{
    tab <- tab.expression()
    if (file.exists(file.path(data.dir,"tables","expression_RNAseq_time_course.tab"))){
      text <- "expression_RNAseq_time_course.tab already exists in the tables folder"
    }
    else {
      write.table(tab,file.path(data.dir,"tables","expression_RNAseq_time_course.tab"),sep="\t",quote=F,row.names=F)
      text <- "Interactions have been saved in expression_RNAseq_time_course.tab"
    }
    
    text
  })
  
  output$write.tab.expression <- renderText({
    write.tab.expression()
  })
  
  
  #################################################################################################
  ## frequence per cluster
  #################################################################################################
  
  freq.TFBS <- reactive({
    message("Preparing TFBS/ChIP frequency per cluster ...")
    homer.tmp <- prepare.shaps()$homer.tmp
    homer.tmp$db.clust <- as.vector(homer.tmp$db.clust)
    homer.tmp$db.clust[!is.na(homer.tmp$cl1.split)] <- homer.tmp$cl1.split[!is.na(homer.tmp$cl1.split)] 
    TFBS <- data.table(homer.tmp[,c(1:416,grep("db.clust",colnames(homer.tmp)))])
    chip <- data.table(homer.tmp[,c(432:466,grep("db.clust",colnames(homer.tmp)))])
    freq.TFBS <- TFBS[,lapply(.SD,function(x) sum(x>0)/length(x)),by=db.clust]
    chip.TFBS <- chip[,lapply(.SD,function(x) sum(x>0)/length(x)),by=db.clust]
    list(TFBS=freq.TFBS,chip=chip.TFBS)
  })
  
  output$heatmap.freq.TFBS <- renderPlot({
    mat <- freq.TFBS()$TFBS
    mat <- as.data.frame(mat)
    row.names(mat)  <- mat$db.clust
    mat <- mat[,-1]
    mat <- mat[,apply(mat,2,function(x){max(x)>0.2})]
    pheatmap(mat,clustering_distance_cols = "binary", clustering_distance_rows = "binary",
             clustering_method = "ward.D",fontsize = 7,color=colorRampPalette(colours()[c(1,142,53,34,35,24)])(100),breaks=seq(0,1,0.01))
  })
  
  output$heatmap.freq.chip <- renderPlot({
    mat <- freq.TFBS()$chip
    mat <- as.data.frame(mat)
    row.names(mat)  <- mat$db.clust
    mat <- mat[,-1]
    mat <- mat[,apply(mat,2,function(x){max(x)>0.1})]
    pheatmap(mat,clustering_distance_cols = "binary", clustering_distance_rows = "binary",
             clustering_method = "ward.D",fontsize = 7,color=colorRampPalette(colours()[c(1,142,53,34,35,24)])(100),breaks=seq(0,1,0.01))
  })
  
  bait.cluster.interaction <- reactive({
    interacting.tSNE <- interacting.tSNE()
    # m[,c("heatmap.info","class","cluster1","ATAC.id1","fragment.id1","ATAC.id2","fragment.id2","cluster2")]
    sense <- interacting.tSNE[grepl("^chr",fragment.id2)]
    antisense <- interacting.tSNE[grepl("^chr",fragment.id1)]
    interacting.tSNE <- rbind(sense[,list(fragment.id1,cluster1,fragment.id2,cluster2,ATAC.id2)],
                              antisense[,list(fragment.id1=fragment.id2,cluster1=cluster2,fragment.id2=fragment.id1,cluster2=cluster1,ATAC.id2=ATAC.id1)])
    interacting.tSNE <- unique(interacting.tSNE)
    if (!"no.ATAC"%in%input$display.no.clust){
      interacting.tSNE <- interacting.tSNE[cluster2!="no.ATAC"]
    }
    if (!"not.in.tSNE"%in%input$display.no.clust){
      interacting.tSNE <- interacting.tSNE[cluster2!="not.in.tSNE"]
    }
    if (!"cluster0"%in%input$display.no.clust){
      interacting.tSNE <- interacting.tSNE[cluster2!="0"]
    }
    clust.assoc <- as.data.frame.matrix(table(interacting.tSNE$fragment.id1,interacting.tSNE$cluster2))
    
    list(clust.assoc=clust.assoc,raw=interacting.tSNE)
  })
  
  output$bait.cluster.interaction.plot <- renderPlot({
    bait.cluster.interaction <- bait.cluster.interaction()$clust.assoc
    if (input$cl.inter=="count"){
      pheatmap(as.data.frame(t(bait.cluster.interaction)),show_rownames = T, 
               show_colnames = F,
               border_color = NA,
               clustering_distance_cols = "binary",
               color = colorRampPalette(c("white","orange","red","brown","black"))(100),
               breaks = seq(0,10,0.1))
    }
    else if (input$cl.inter=="frequency"){
      tot <- apply(bait.cluster.interaction,1,sum)
      bait.cluster.interaction <- bait.cluster.interaction/tot
      pheatmap(as.data.frame(t(bait.cluster.interaction)),show_rownames = T, 
               show_colnames = F,
               border_color = NA,
               clustering_distance_cols = "binary",
               color = colorRampPalette(c("white","orange","red","brown","black"))(100),
               breaks = seq(0,1,0.01))
    }
  })
  
  output$choose.cls <- renderUI({
    mat <- bait.cluster.interaction()$clust.assoc
    checkboxGroupInput("choose.cl","Select clusters having other fragments",choices=c("all",colnames(mat)),inline = T)
  })
  
  show.bait.w.other.cl <- eventReactive(input$go.show.me.genes,{
    mat <- bait.cluster.interaction()$clust.assoc
    print(dim(mat))
    cl <- input$choose.cl
    print(cl)
    if (length(cl)>1){
      select.rows <- apply(mat[,cl],1,function(x,opt){
        if (opt=="AND"){
          f <- sum(x>0) == length(x)
        }
        else {
          f <- sum(x>0) > 0
        }
        
      },input$AND.OR.cl.choice)
    }
    else if (length(cl)==1 & cl!="all"){
      select.rows <- mat[,cl]>0
    }
    else if (length(cl)==1 & cl=="all"){
      select.rows <- rep(T,ncol(mat))
    }
    else {
      select.rows <- F
    }
    text.show <- paste("The cluster(s)",paste(cl,collapse = ", "), "contain(s) fragments interacting with", sum(select.rows), "genes: ",paste(row.names(mat)[select.rows],collapse=", "))
    if (input$cl.inter=="frequency"){
      tot <- apply(mat,1,sum)
      mat <- mat/tot
    }
    if (cl[1]!="all") {
      mat <- mat[select.rows,cl,drop=F]
    }
    mat$SYMBOL <- row.names(mat)
    mat <- mat[,sort(colnames(mat),decreasing = T)]
    print(select.rows)
    list(text.show=text.show,mat=mat)
  })
  
  output$show.bait.w.other.cl <- renderText({
    select.rows <- show.bait.w.other.cl()$text.show
  })
  
  output$show.bait.w.other.matrix <- renderTable({
    select.rows <- show.bait.w.other.cl()$mat
  })
  
  tab.cl.frag.int <- eventReactive(input$download.sel.tab,{
    select.rows <- show.bait.w.other.cl()$mat
    cl <- input$choose.cl
    write.table(select.rows,file.path(data.dir,"tables",paste0("selected_genes_with_int_in_cl_",paste(cl,collapse="_"),".tab")),sep="\t",quote=F)
    text <- paste0("File ",paste0("selected_genes_with_int_in_cl_",paste(cl,collapse="_"),".tab")," has been saved.")
  })  
  
  output$tab.cl.frag.int <- renderText({
    tab.cl.frag.int()
  })
  
  #################################################################################################
  ## CHiC plots 
  #################################################################################################
  
  int.fragments <- reactive({
    int.fragments <- all.interactions[,list(id1,id2,name1,name2,type1,type2,type,diff.int=extend.classif,TAD,distance,direction)]
    
    m <- match(int.fragments$name1,RNAseq$ID)
    int.fragments$expr1 <- RNAseq$status[m]
    int.fragments$logFC1 <- RNAseq$logFC[m]
    
    m <- match(int.fragments$name2,RNAseq$ID)
    int.fragments$expr2 <- RNAseq$status[m]
    int.fragments$logFC2 <- RNAseq$logFC[m]
    
    #int.fragments <- int.fragments[!is.na(expr1)]
    #int.fragments <- int.fragments[type=="bo" | (type=="bb" & !is.na(expr2))]
    
    int.fragments <- rbind(int.fragments[,list(id=id1,name=name1,type=type1,B2B=type,diff.int,TAD,distance,direction)],int.fragments[,list(id=id2,name=name2,type=type2,B2B=type,diff.int,TAD,distance,direction=-direction)])
  })
  
  
  ranges2 <- reactiveValues(x = NULL, y = NULL)
  
  observe({
    brush <- input$plot2_zoom
    if (!is.null(brush)) {
      ranges2$x <- c(brush$xmin, brush$xmax)
      ranges2$y <- c(brush$ymin, brush$ymax)
      
    } else {
      ranges2$x <- NULL
      ranges2$y <- NULL
    }
  })
  
  observeEvent(input$zoomin,{
    if (!is.null(input$plot2_zoom)) {
      ranges2$x <- sort(c(ranges2$x[1]+(as.numeric(input$kb.zoom)*1000),ranges2$x[2]-(as.numeric(input$kb.zoom)*1000)))
      ranges2$y <- ranges2$y
    } else {
      ranges2$x <- NULL
      ranges2$y <- NULL
    }
    
  })
  observeEvent(input$zoomout,{
    if (!is.null(input$plot2_zoom)) {
      ranges2$x <- sort(c(ranges2$x[1]-(as.numeric(input$kb.zoom)*1000),ranges2$x[2]+(as.numeric(input$kb.zoom)*1000)))
      ranges2$y <- ranges2$y
    } else {
      ranges2$x <- NULL
      ranges2$y <- NULL
    }
  })
  
  observeEvent(input$Go.left,{
    if (!is.null(input$plot2_zoom)) {
      ranges2$x <- sort(c(ranges2$x[1]-(as.numeric(input$kb.move)*1000),ranges2$x[2]-(as.numeric(input$kb.move)*1000)))
      ranges2$y <- ranges2$y
    } else {
      ranges2$x <- NULL
      ranges2$y <- NULL
    }
  })
  
  observeEvent(input$Go.right,{
    if (!is.null(input$plot2_zoom)) {
      ranges2$x <- sort(c(ranges2$x[1]+(as.numeric(input$kb.move)*1000),ranges2$x[2]+(as.numeric(input$kb.move)*1000)))
      ranges2$y <- ranges2$y
    } else {
      ranges2$x <- NULL
      ranges2$y <- NULL
    }
  })
  
  check.region.of <- eventReactive(input$ckeck.exists,{
    coords <- input$coords.TADs
    if (grepl("^chr",coords)){
      chrom <- unlist(sapply(strsplit(coords,":"),"[[",1))
      coords <- unlist(sapply(strsplit(coords,":"),"[[",2))
      s <- as.numeric(unlist(sapply(strsplit(coords,"-"),"[[",1)))
      e <- as.numeric(unlist(sapply(strsplit(coords,"-"),"[[",2)))
      inter <- all.interactions[chr1==chrom]
      inter <- inter[(locus1>=s & locus1<=e) | (locus2>=s & locus2<=e)]
      size.area <- round(abs(min(c(inter$locus1,inter$locus2))-max(c(inter$locus1,inter$locus2)))/10^6,2)
      area1 <- paste0(chrom,":",min(c(inter$locus1,inter$locus2)),"-",max(c(inter$locus1,inter$locus2)))
      text <- paste("The area (",area1,") is around",size.area,"Mb and contains",nrow(inter),"stable interactions.")
      if (sum(inter$extend.classif!="stable")>0){
        inter2 <- all.interactions[chr1==chrom & extend.classif!="stable"]
        inter2 <- inter2[(locus1>=s & locus1<=e) | (locus2>=s & locus2<=e)]
        size.area2 <- round(abs(min(c(inter2$locus1,inter2$locus2))-max(c(inter2$locus1,inter2$locus2)))/10^6,2)
        area2 <- paste0(chrom,":",min(c(inter2$locus1,inter2$locus2)),"-",max(c(inter2$locus1,inter2$locus2)))
        
        text <- paste("The area (",area1,") can be up to",size.area,"Mb if displaying all interactions (",nrow(inter),") or to",size.area2,"Mb (",area2,") if selecting differential (",sum(inter$extend.classif!="stable"),").")
        
      }
    }
    else {
      SYMBOL <- coords
      if (SYMBOL%in%c(all.interactions$name1,all.interactions$name2)){
        
        inter <- all.interactions[name1==SYMBOL | name2==SYMBOL]
        chrom <- unique(inter$chr1)
        size.area <- round(abs(min(c(inter$locus1,inter$locus2))-max(c(inter$locus1,inter$locus2)))/10^6,2)
        area1 <- paste0(chrom,":",min(c(inter$locus1,inter$locus2)),"-",max(c(inter$locus1,inter$locus2)))
        text <- paste("The area (",area1,") is around",size.area,"Mb and contains",nrow(inter),"stable interactions.")
        if (sum(inter$extend.classif!="stable")>0){
          inter2 <- all.interactions[(name1==SYMBOL | name2==SYMBOL) & extend.classif!="stable"]
          size.area2 <- round(abs(min(c(inter2$locus1,inter2$locus2))-max(c(inter2$locus1,inter2$locus2)))/10^6,2)
          area2 <- paste0(chrom,":",min(c(inter2$locus1,inter2$locus2)),"-",max(c(inter2$locus1,inter2$locus2)))
          text <- paste("The area (",area1,") can be up to",size.area,"Mb if displaying all interactions (",nrow(inter),") or to",size.area2,"Mb (",area2,") if selecting differential (",sum(inter$extend.classif!="stable"),").")
        }
      }
      else {
        text <- paste(SYMBOL, "does not have any related data. Applying would return an error, try another identifier.")
      }
    }
    text
  })
  
  output$region.exists <- renderText({
    check.region.of()
  })
  
  plot.TAD.summ <- eventReactive(input$go.CHiC.plot.TAD,{
    message("Preparing interaction plot ...")
    coords <- input$coords.TADs
    message(coords,"taken into account ...")
    if (grepl("^chr",coords)){
      chrom <- unlist(sapply(strsplit(coords,":"),"[[",1))
      coords <- unlist(sapply(strsplit(coords,":"),"[[",2))
      s <- as.numeric(unlist(sapply(strsplit(coords,"-"),"[[",1)))
      e <- as.numeric(unlist(sapply(strsplit(coords,"-"),"[[",2)))
      inter <- all.interactions[chr1==chrom]
      inter <- inter[(locus1>=s & locus1<=e) | (locus2>=s & locus2<=e)]
    }
    else {
      SYMBOL <- coords
      selected.interactions <- all.interactions[name1==SYMBOL | name2==SYMBOL]
      chrom <- as.vector(seqnames(all.fragments))[all.fragments$gene==SYMBOL | all.fragments$IDs==SYMBOL]
      start.frag <- as.vector(start(all.fragments))[all.fragments$gene==SYMBOL | all.fragments$IDs==SYMBOL]
      s <- min(c(selected.interactions$locus1,selected.interactions$locus2))
      e <- max(c(selected.interactions$locus1,selected.interactions$locus2))
      inter <- selected.interactions
    }
    inter$end1 <- end(all.fragments)[match(inter$id1,all.fragments$IDs)]
    inter$end2 <- end(all.fragments)[match(inter$id2,all.fragments$IDs)]
    int.fragments.stats <- per.fragment.interactions.total()
    
    to.plot <- int.fragments.stats[total>1]
    if (input$arc.plot=="all"){
      s <- min(c(inter$locus1,inter$locus2))
      e <- max(c(inter$end1,inter$end2))
      g <- ggplot(to.plot[chr==chrom],aes(x=start,y=1)) + 
        geom_segment(data=as.data.frame(tads_borders[as.vector(seqnames(tads_borders))=="chr1"]),aes(x=start,xend=start,y=-1,yend=1),linetype="dashed",color=colours()[613]) + 
        geom_segment(aes(xend=start,yend=0.75),size=0.5,colour=ifelse(to.plot[chr==chrom]$type=="o","darkgrey","black")) + 
        xlim(s,e) + ylim(-1,1.5) + theme_bw()
      coords <- as.data.frame(t(apply(inter[,list(locus1,locus2)],1,sort)))
      inter$locus1 <- coords$V1
      inter$locus2 <- coords$V2
      g <- g + geom_curve(curvature = 1,lineend = 'butt',data=inter,aes(x=locus1,xend=locus2,y=0,yend=0),linetype="dashed",size=0.2,color="grey",alpha=1)#ifelse(inter$extend.classif=="gain","red",colours()[124]))
      if(sum(inter$extend.classif!="stable")>0){
        g <- g + geom_curve(curvature = 1,lineend = 'butt',data=inter[extend.classif!="stable"],aes(x=locus1,xend=locus2,y=0,yend=0),linetype="dashed",color=ifelse(inter[extend.classif!="stable"]$extend.classif=="gain","red",colours()[124]))
      }
    }
    if(input$arc.plot=="differential" & sum(inter$extend.classif!="stable")>0){
      inter <- inter[extend.classif!="stable"]
      s <- min(c(inter$locus1,inter$locus2))
      e <- max(c(inter$end1,inter$end2))
      g <- ggplot(to.plot[chr==chrom],aes(x=start,y=1)) + 
        geom_segment(data=as.data.frame(tads_borders[as.vector(seqnames(tads_borders))=="chr1"]),aes(x=start,xend=start,y=-1,yend=1),linetype="dashed",size=1.2,color=colours()[613]) + 
        geom_segment(aes(xend=start,yend=0.75),size=0.5,colour=ifelse(to.plot[chr==chrom]$type=="o","darkgrey","black")) + 
        theme_bw()
      coords <- as.data.frame(t(apply(inter[,list(locus1,locus2)],1,sort)))
      inter$locus1 <- coords[,1]
      inter$locus2 <- coords[,2]
      g <- g + geom_curve(curvature = 1,lineend = 'butt',data=inter[extend.classif!="stable"],aes(x=locus1,xend=locus2,y=0,yend=0),linetype="dashed",color=ifelse(inter[extend.classif!="stable"]$extend.classif=="gain","red",colours()[124]))
    }
    
    else if(input$arc.plot=="differential" & sum(inter$extend.classif!="stable")==0){
      s <- min(c(inter$locus1,inter$locus2))
      e <- max(c(inter$end1,inter$end2))
      g <- ggplot(to.plot[chr==chrom],aes(x=start,y=1)) + 
        geom_segment(data=as.data.frame(tads_borders[as.vector(seqnames(tads_borders))=="chr1"]),aes(x=start,xend=start,y=-1,yend=1),linetype="dashed",color=colours()[613]) + 
        geom_segment(aes(xend=start,yend=0.75),size=0.5,colour=ifelse(to.plot[chr==chrom]$type=="o","darkgrey","black")) + 
        xlim(s,e) + ylim(-1,1.5) + theme_bw()
      coords <- as.data.frame(t(apply(inter[,list(locus1,locus2)],1,sort)))
      inter$locus1 <- coords$V1
      inter$locus2 <- coords$V2
      g <- g + ggtitle("No differential interaction with this fragment or in this region. Showing stable ones.") + geom_curve(curvature = 1,lineend = 'butt',data=inter,aes(x=locus1,xend=locus2,y=0,yend=0),linetype="dashed",size=0.2,color="grey",alpha=1)#ifelse(inter$extend.classif=="gain","red",colours()[124]))
    }
    
    symbols.text <- int.fragments.stats[chr==chrom & start>=s & start<=e & type=="b"]
    
    g <- g + scale_colour_manual(values=c(Th1=colours()[54],Th2=colours()[613],up="red",down=colours()[131],stable="black",no.data="darkgrey",opening="red",closing=colours()[124],G0=colours()[124],G1="red")) +
      scale_fill_manual(values=c(up="red",down=colours()[131],stable="black",no.data="darkgrey",opening="red",closing=colours()[124]))
    
    genes <- data.table(as.data.frame(refSeq[as.vector(seqnames(refSeq))==chrom & ((start(refSeq)>s & start(refSeq)<e) | (end(refSeq)>s & end(refSeq)<e))]))
    if (nrow(genes)>0){
      s <- min(c(min(c(genes$start,genes$end)),s)) - 100
      e <- max(c(max(c(genes$start,genes$end)),e)) + 100
      
      genes$status <- RNAseq[match(genes$SYMBOL,ID)]$status
      genes$status[is.na(genes$status)] <- "no.data"
      genes$directionality <- 0
      g <- g + geom_rect(data=genes,aes(xmin=start,xmax=end,ymin=-0.05,ymax=0.05,fill=status,colour=status),alpha=0.5)
      
    }
    if (nrow(symbols.text)>0){
      symbols.text[total<=10]$directionality <- 0
      symbols.text$status <- RNAseq[match(symbols.text$name,ID)]$status
      symbols.text$status[is.na(symbols.text$status)] <- "no.data"
      g <- g + geom_text_repel(data=symbols.text,aes(x=start,y=0.75,label=name,color=status),angle=45) 
    }
    ATAC.select <- as.data.frame(peaks.ATAC)
    ATAC.select <- ATAC.select[ATAC.select$seqnames==chrom & ATAC.select$start>s-100 & ATAC.select$end<e+100,]
    ATAC.select$center <- (ATAC.select$start + ATAC.select$end)/2
    frags <- as.data.frame(all.fragments[all.fragments$IDs%in%c(inter$id1,inter$id2)])
    g <- g + geom_hline(yintercept = 0,size=0.3,linetype="dashed") + 
      geom_segment(data=frags,aes(x=start+10,xend=end-10,y=1,yend=1),alpha=0.8,size=5) + 
      geom_label(aes(x=s,y=1,label="fragments")) + geom_label(aes(x=s,y=1.2,label="ATAC"))
    
    y.pos <- 1.35
    
    if (nrow(ATAC.select)>0){
      g <- g + geom_segment(data=ATAC.select,aes(x=center,xend=center,y=0,yend=1.2,colour=status),size=0.1,linetype="dashed") + geom_segment(data=ATAC.select,aes(x=start,xend=end,y=1.2,yend=1.2,colour=status),size=5) 
      g2 <- ggplot(ATAC.select,aes(x=start,xend=end,y=1.2,yend=1.2)) + geom_segment(aes(colour=status),size=5)
    }
    else {
      ATAC.select <- data.table(start=1,end=2,status="stable")
      g2 <- ggplot(ATAC.select,aes(x=start,xend=end,y=1.2,yend=1.2)) + geom_segment(aes(colour=status),size=5)
    }
    
    format.mark <- function(H3,m,s,e,chrom){
      mark <- H3[H3$mark==m]
      mark <- data.table(as.data.frame(mark[as.vector(seqnames(mark))==chrom & ((start(mark)>s & start(mark)<e) | (end(mark)>s & end(mark)<e))]))
    }
    
    yg2 <- "ATAC"
    ys <- 1.2
    YY1 <- data.table(as.data.frame(YY1[as.vector(seqnames(YY1))==chrom & ((start(YY1)>s & start(YY1)<e) | (end(YY1)>s & end(YY1)<e))]))
    
    if (nrow(YY1)>0){
      yg2 <- c(yg2,"YY1")
      y.tmpYY1 <- y.pos
      ys <- c(ys,y.pos)
      s <- min(c(min(c(YY1$start,YY1$end)),s)) - 100
      e <- max(c(max(c(YY1$start,YY1$end)),e)) + 100
      
      g2 <- g2 + geom_segment(data=YY1,aes(x=start,xend=end,y=y.tmpYY1,yend=y.tmpYY1,color=diff),alpha=0.5,size=5) + geom_label(aes(x=s,y=y.tmpYY1,label="YY1"))
      y.pos <- y.pos + 0.15
    }
    
    
    cFos <- reduce(GSE116695.peaks$cFos_3)
    cFos <- data.table(as.data.frame(cFos[as.vector(seqnames(cFos))==chrom & ((start(cFos)>s & start(cFos)<e) | (end(cFos)>s & end(cFos)<e))]))
    
    if (nrow(cFos)>0){
      yg2 <- c(yg2,"cFos")
      ys <- c(ys,y.pos)
      
      y.tmpFOS <- y.pos
      s <- min(c(min(c(cFos$start,cFos$end)),s)) - 100
      e <- max(c(max(c(cFos$start,cFos$end)),e)) + 100
      
      g2 <- g2 + geom_segment(data=cFos,aes(x=start,xend=end,y=y.tmpFOS,yend=y.tmpFOS),alpha=0.5,size=5,color=colours()[465]) + geom_label(aes(x=s,y=y.tmpFOS,label="cFos"))
      y.pos <- y.pos + 0.15
    }
    
    JunB <- reduce(GSE116695.peaks$JunB_2)
    JunB <- data.table(as.data.frame(JunB[as.vector(seqnames(JunB))==chrom & ((start(JunB)>s & start(JunB)<e) | (end(JunB)>s & end(JunB)<e))]))
    
    if (nrow(JunB)>0){
      yg2 <- c(yg2,"JunB")
      ys <- c(ys,y.pos)
      
      y.tmpJun <- y.pos
      s <- min(c(min(c(JunB$start,JunB$end)),s)) - 100
      e <- max(c(max(c(JunB$start,JunB$end)),e)) + 100
      
      g2 <- g2 + geom_segment(data=JunB,aes(x=start,xend=end,y=y.tmpJun,yend=y.tmpJun),alpha=0.5,size=5,color=colours()[465]) + geom_label(aes(x=s,y=y.tmpJun,label="JunB"))
      y.pos <- y.pos + 0.15
    }
    
    cMyc <- reduce(GSE116695.peaks$cMyc_2)
    cMyc <- data.table(as.data.frame(cMyc[as.vector(seqnames(cMyc))==chrom & ((start(cMyc)>s & start(cMyc)<e) | (end(cMyc)>s & end(cMyc)<e))]))
    
    if (nrow(cMyc)>0){
      yg2 <- c(yg2,"cMyc")
      ys <- c(ys,y.pos)
      
      y.tmpMYC <- y.pos
      s <- min(c(min(c(cMyc$start,cMyc$end)),s)) - 100
      e <- max(c(max(c(cMyc$start,cMyc$end)),e)) + 100
      
      g2 <- g2 + geom_segment(data=cMyc,aes(x=start,xend=end,y=y.tmpMYC,yend=y.tmpMYC),alpha=0.5,size=5,color=colours()[465]) + geom_label(aes(x=s,y=y.tmpMYC,label="cMyc"))
      y.pos <- y.pos + 0.15
    }
    
    p50 <- reduce(c(GSE116695.peaks$p50_1,GSE116695.peaks$p50_3))
    p50 <- data.table(as.data.frame(p50[as.vector(seqnames(p50))==chrom & ((start(p50)>s & start(p50)<e) | (end(p50)>s & end(p50)<e))]))
    
    if (nrow(p50)>0){
      yg2 <- c(yg2,"p50")
      ys <- c(ys,y.pos)
      
      y.tmpp50 <- y.pos
      s <- min(c(min(c(p50$start,p50$end)),s)) - 100
      e <- max(c(max(c(p50$start,p50$end)),e)) + 100
      
      g2 <- g2 + geom_segment(data=p50,aes(x=start,xend=end,y=y.tmpp50,yend=y.tmpp50),alpha=0.5,size=5,color=colours()[465]) + geom_label(aes(x=s,y=y.tmpp50,label="p50"))
      y.pos <- y.pos + 0.15
    }
    
    NFAT1 <- reduce(c(GSE116695.peaks$NFAT1_1,GSE116695.peaks$NFAT1_3))
    NFAT1 <- data.table(as.data.frame(NFAT1[as.vector(seqnames(NFAT1))==chrom & ((start(NFAT1)>s & start(NFAT1)<e) | (end(NFAT1)>s & end(NFAT1)<e))]))
    
    if (nrow(NFAT1)>0){
      yg2 <- c(yg2,"NFAT1")
      ys <- c(ys,y.pos)
      
      y.tmpNFAT1 <- y.pos
      s <- min(c(min(c(NFAT1$start,NFAT1$end)),s)) - 100
      e <- max(c(max(c(NFAT1$start,NFAT1$end)),e)) + 100
      
      g2 <- g2 + geom_segment(data=NFAT1,aes(x=start,xend=end,y=y.tmpNFAT1,yend=y.tmpNFAT1),alpha=0.5,size=5,color=colours()[465]) + geom_label(aes(x=s,y=y.tmpNFAT1,label="NFAT1"))
      y.pos <- y.pos + 0.15
    }
    
    NFAT2 <- reduce(c(GSE116695.peaks$NFAT2_1,GSE116695.peaks$NFAT2_3))
    NFAT2 <- data.table(as.data.frame(NFAT2[as.vector(seqnames(NFAT2))==chrom & ((start(NFAT2)>s & start(NFAT2)<e) | (end(NFAT2)>s & end(NFAT2)<e))]))
    
    if (nrow(NFAT2)>0){
      yg2 <- c(yg2,"NFAT2")
      ys <- c(ys,y.pos)
      
      y.tmpNFAT2 <- y.pos
      s <- min(c(min(c(NFAT2$start,NFAT2$end)),s)) - 100
      e <- max(c(max(c(NFAT2$start,NFAT2$end)),e)) + 100
      
      g2 <- g2 + geom_segment(data=NFAT2,aes(x=start,xend=end,y=y.tmpNFAT2,yend=y.tmpNFAT2),alpha=0.5,size=5,color=colours()[465]) + geom_label(aes(x=s,y=y.tmpNFAT2,label="NFAT2"))
      y.pos <- y.pos + 0.15
    }
    
    RApteb_Th1 <- reduce(GSE62486.peaks$Th1_reactivated_pTEFb)
    RApteb_Th1$cond <- "Th1"
    RApteb_Th2 <- reduce(GSE62486.peaks$Th2_reactivated_pTEFb)
    RApteb_Th2$cond <- "Th2"
    RApteb <- c(RApteb_Th1,RApteb_Th2)
    RApteb <- data.table(as.data.frame(RApteb[as.vector(seqnames(RApteb))==chrom & ((start(RApteb)>s & start(RApteb)<e) | (end(RApteb)>s & end(RApteb)<e))]))
    
    if (nrow(RApteb)>0){
      yg2 <- c(yg2,"RA_pTEFb")
      ys <- c(ys,y.pos)
      
      y.tmpRAptefb <- y.pos
      s <- min(c(min(c(RApteb$start,RApteb$end)),s)) - 100
      e <- max(c(max(c(RApteb$start,RApteb$end)),e)) + 100
      
      g2 <- g2 + geom_segment(data=RApteb,aes(x=start,xend=end,y=y.tmpRAptefb,yend=y.tmpRAptefb,color=cond),alpha=0.5,size=5) + geom_label(aes(x=s,y=y.tmpRAptefb,label="RA_pTEFb"))
      y.pos <- y.pos + 0.15
    }
    
    USpteb_Th1 <- reduce(GSE62486.peaks$`Th1_US_PTEF-b`)
    USpteb_Th1$cond <- "Th1"
    USpteb_Th2 <- reduce(GSE62486.peaks$`Th2_US_PTEF-b`)
    USpteb_Th2$cond <- "Th2"
    USpteb <- c(USpteb_Th1,USpteb_Th2)
    USpteb <- data.table(as.data.frame(USpteb[as.vector(seqnames(USpteb))==chrom & ((start(USpteb)>s & start(USpteb)<e) | (end(USpteb)>s & end(USpteb)<e))]))
    
    if (nrow(USpteb)>0){
      yg2 <- c(yg2,"US_pTEFb")
      ys <- c(ys,y.pos)
      
      y.tmpUSptefb <- y.pos
      s <- min(c(min(c(USpteb$start,USpteb$end)),s)) - 100
      e <- max(c(max(c(USpteb$start,USpteb$end)),e)) + 100
      
      g2 <- g2 + geom_segment(data=USpteb,aes(x=start,xend=end,y=y.tmpUSptefb,yend=y.tmpUSptefb,color=cond),alpha=0.5,size=5) + geom_label(aes(x=s,y=y.tmpUSptefb,label="US_pTEFb"))
      y.pos <- y.pos + 0.15
    }
    
    USpteb_DMSO <- reduce(GSE62486.peaks$`Th1_DMSO_P-TEFb`)
    USpteb_DMSO$cond <- "DMSO"
    USpteb_BAY <- reduce(GSE62486.peaks$`Th1_BAY_P-TEFb`)
    USpteb_BAY$cond <- "BAY"
    pteb <- c(USpteb_DMSO,USpteb_BAY)
    pteb <- data.table(as.data.frame(pteb[as.vector(seqnames(pteb))==chrom & ((start(pteb)>s & start(pteb)<e) | (end(pteb)>s & end(pteb)<e))]))
    
    if (nrow(pteb)>0){
      yg2 <- c(yg2,"pTEFb")
      ys <- c(ys,y.pos)
      
      y.tmpptefb <- y.pos
      s <- min(c(min(c(pteb$start,pteb$end)),s)) - 100
      e <- max(c(max(c(pteb$start,pteb$end)),e)) + 100
      
      g2 <- g2 + geom_segment(data=pteb,aes(x=start,xend=end,y=y.tmpptefb,yend=y.tmpptefb,color=cond),alpha=0.5,size=5) + geom_label(aes(x=s,y=y.tmpptefb,label="pTEFb"))
      y.pos <- y.pos + 0.15
    }
    
    USmed1_DMSO <- reduce(GSE62486.peaks$Th1_RA_Med1)
    USmed1_DMSO$cond <- "RA"
    USmed1_BAY <- reduce(GSE62486.peaks$Th1_BAY_Med1)
    USmed1_BAY$cond <- "BAY"
    med1 <- c(USmed1_DMSO,USmed1_BAY)
    med1 <- data.table(as.data.frame(med1[as.vector(seqnames(med1))==chrom & ((start(med1)>s & start(med1)<e) | (end(med1)>s & end(med1)<e))]))
    
    if (nrow(med1)>0){
      yg2 <- c(yg2,"Med1")
      ys <- c(ys,y.pos)
      
      y.tmpmed1 <- y.pos
      s <- min(c(min(c(med1$start,med1$end)),s)) - 100
      e <- max(c(max(c(med1$start,med1$end)),e)) + 100
      
      g2 <- g2 + geom_segment(data=med1,aes(x=start,xend=end,y=y.tmpmed1,yend=y.tmpmed1,color=cond),alpha=0.5,size=5) + geom_label(aes(x=s,y=y.tmpmed1,label="Med1"))
      y.pos <- y.pos + 0.15
    }
    
    tbet<- reduce(GSE62486.peaks$`Th1_WT_T-bet`)
    
    tbet <- data.table(as.data.frame(tbet[as.vector(seqnames(tbet))==chrom & ((start(tbet)>s & start(tbet)<e) | (end(tbet)>s & end(tbet)<e))]))
    
    if (nrow(tbet)>0){
      yg2 <- c(yg2,"T-bet")
      ys <- c(ys,y.pos)
      
      y.tmptbet <- y.pos
      s <- min(c(min(c(tbet$start,tbet$end)),s)) - 100
      e <- max(c(max(c(tbet$start,tbet$end)),e)) + 100
      
      g2 <- g2 + geom_segment(data=tbet,aes(x=start,xend=end,y=y.tmptbet,yend=y.tmptbet),alpha=0.5,size=5) + geom_label(aes(x=s,y=y.tmptbet,label="T-bet"))
      y.pos <- y.pos + 0.15
    }
    
    USAff4_DMSO <- reduce(GSE62486.peaks$Th1_RA_Aff4)
    USAff4_DMSO$cond <- "RA"
    USAff4_BAY <- reduce(GSE62486.peaks$Th1_BAY_Aff4)
    USAff4_BAY$cond <- "BAY"
    Aff4 <- c(USAff4_DMSO,USAff4_BAY)
    Aff4 <- data.table(as.data.frame(Aff4[as.vector(seqnames(Aff4))==chrom & ((start(Aff4)>s & start(Aff4)<e) | (end(Aff4)>s & end(Aff4)<e))]))
    
    if (nrow(Aff4)>0){
      yg2 <- c(yg2,"Aff4")
      ys <- c(ys,y.pos)
      
      y.tmpAff4 <- y.pos
      s <- min(c(min(c(Aff4$start,Aff4$end)),s)) - 100
      e <- max(c(max(c(Aff4$start,Aff4$end)),e)) + 100
      
      g2 <- g2 + geom_segment(data=Aff4,aes(x=start,xend=end,y=y.tmpAff4,yend=y.tmpAff4,color=cond),alpha=0.5,size=5) + geom_label(aes(x=s,y=y.tmpAff4,label="Aff4"))
      y.pos <- y.pos + 0.15
    }
    
    USp65_DMSO <- reduce(GSE62486.peaks$Th1_RA_p65)
    USp65_DMSO$cond <- "RA"
    USp65_BAY <- reduce(GSE62486.peaks$Th1_BAY_p65)
    USp65_BAY$cond <- "BAY"
    p65 <- c(USp65_DMSO,USp65_BAY)
    p65 <- data.table(as.data.frame(p65[as.vector(seqnames(p65))==chrom & ((start(p65)>s & start(p65)<e) | (end(p65)>s & end(p65)<e))]))
    
    if (nrow(p65)>0){
      yg2 <- c(yg2,"p65")
      ys <- c(ys,y.pos)
      
      y.tmpp65 <- y.pos
      s <- min(c(min(c(p65$start,p65$end)),s)) - 100
      e <- max(c(max(c(p65$start,p65$end)),e)) + 100
      
      g2 <- g2 + geom_segment(data=p65,aes(x=start,xend=end,y=y.tmpp65,yend=y.tmpp65,color=cond),alpha=0.5,size=5) + geom_label(aes(x=s,y=y.tmpp65,label="p65"))
      y.pos <- y.pos + 0.15
    }
    
    USBrd4_DMSO <- reduce(GSE62486.peaks$Th1_RA_Brd4)
    USBrd4_DMSO$cond <- "RA"
    USBrd4_BAY <- reduce(GSE62486.peaks$Th1_BAY_Brd4)
    USBrd4_BAY$cond <- "BAY"
    Brd4 <- c(USBrd4_DMSO,USBrd4_BAY)
    Brd4 <- data.table(as.data.frame(Brd4[as.vector(seqnames(Brd4))==chrom & ((start(Brd4)>s & start(Brd4)<e) | (end(Brd4)>s & end(Brd4)<e))]))
    
    if (nrow(Brd4)>0){
      yg2 <- c(yg2,"Brd4")
      ys <- c(ys,y.pos)
      
      y.tmpBrd4 <- y.pos
      s <- min(c(min(c(Brd4$start,Brd4$end)),s)) - 100
      e <- max(c(max(c(Brd4$start,Brd4$end)),e)) + 100
      
      g2 <- g2 + geom_segment(data=Brd4,aes(x=start,xend=end,y=y.tmpBrd4,yend=y.tmpBrd4,color=cond),alpha=0.5,size=5) + geom_label(aes(x=s,y=y.tmpBrd4,label="Brd4"))
      y.pos <- y.pos + 0.15
    }
    
    m <- "H3K27Ac"; mark <-   format.mark(H3,m,s,e,chrom)
    if (nrow(mark)>0){
      yg2 <- c(yg2,"H3K27Ac")
      ys <- c(ys,y.pos)
      
      y.tmp1 <- y.pos
      m.temp1 <- m
      mark.tmp1 <- mark
      s <- min(c(min(c(mark.tmp1$start,mark.tmp1$end)),s)) - 100
      e <- max(c(max(c(mark.tmp1$start,mark.tmp1$end)),e)) + 100
      
      g2 <- g2 + geom_segment(data=mark.tmp1,aes(x=start,xend=end,y=y.tmp1,yend=y.tmp1,colour=cond),alpha=0.5,size=5) + geom_label(aes(x=s,y=y.tmp1,label=m.temp1))
      y.pos <- y.pos + 0.15
    }
    
    k27ac <- reduce(GSE116695.peaks$K27ac_1)
    k27ac$cond <- "5h"
    k27ac_naive <- reduce(GSE116695.peaks$Naive_H3K27Ac)
    k27ac_naive$cond <- "0h"
    k27ac_24 <- reduce(GSE116695.peaks$E24h_K27ac)
    k27ac_24$cond <- "24h"
    k27ac <- c(k27ac,k27ac_naive,k27ac_24)
    k27ac <- data.table(as.data.frame(k27ac[as.vector(seqnames(k27ac))==chrom & ((start(k27ac)>s & start(k27ac)<e) | (end(k27ac)>s & end(k27ac)<e))]))
    
    if (nrow(k27ac)>0){
      yg2 <- c(yg2,"k27ac")
      ys <- c(ys,y.pos)
      
      y.tmp27 <- y.pos
      s <- min(c(min(c(k27ac$start,k27ac$end)),s)) - 100
      e <- max(c(max(c(k27ac$start,k27ac$end)),e)) + 100
      
      g2 <- g2 + geom_segment(data=k27ac,aes(x=start,xend=end,y=y.tmp27,yend=y.tmp27,color=cond),alpha=0.5,size=5) + geom_label(aes(x=s,y=y.tmp27,label="k27ac"))
      y.pos <- y.pos + 0.15
    }
    
    
    k27ac_Th1 <- reduce(GSE62486.peaks$Th1_PI_K27ac)
    k27ac_Th1$cond <- "Th1"
    k27ac_Th2 <- reduce(GSE62486.peaks$Th2_PI_K27ac)
    k27ac_Th2$cond <- "Th2"
    PI_k27ac <- c(k27ac_Th1,k27ac_Th2)
    PI_k27ac <- data.table(as.data.frame(PI_k27ac[as.vector(seqnames(PI_k27ac))==chrom & ((start(PI_k27ac)>s & start(PI_k27ac)<e) | (end(PI_k27ac)>s & end(PI_k27ac)<e))]))
    
    if (nrow(PI_k27ac)>0){
      yg2 <- c(yg2,"PI_k27ac")
      ys <- c(ys,y.pos)
      
      y.tmp27PI <- y.pos
      s <- min(c(min(c(PI_k27ac$start,PI_k27ac$end)),s)) - 100
      e <- max(c(max(c(PI_k27ac$start,PI_k27ac$end)),e)) + 100
      
      g2 <- g2 + geom_segment(data=PI_k27ac,aes(x=start,xend=end,y=y.tmp27PI,yend=y.tmp27PI,color=cond),alpha=0.5,size=5) + geom_label(aes(x=s,y=y.tmp27PI,label="PI_k27ac"))
      y.pos <- y.pos + 0.15
    }
    
    
    m <- "H3K4me1"; mark <-   format.mark(H3,m,s,e,chrom)
    if (nrow(mark)>0){
      yg2 <- c(yg2,"H3K4me1")
      ys <- c(ys,y.pos)
      
      y.tmp2 <- y.pos
      m.temp2 <- m
      mark.tmp2 <- mark
      s <- min(c(min(c(mark.tmp2$start,mark.tmp2$end)),s)) - 100
      e <- max(c(max(c(mark.tmp2$start,mark.tmp2$end)),e)) + 100
      
      g2 <- g2 + geom_segment(data=mark.tmp2,aes(x=start,xend=end,y=y.tmp2,yend=y.tmp2,colour=cond),alpha=0.5,size=5) + geom_label(aes(x=s,y=y.tmp2,label=m.temp2))
      y.pos <- y.pos + 0.15
    }
    m <- "H3K4me3"; mark <-   format.mark(H3,m,s,e,chrom)
    if (nrow(mark)>0){
      yg2 <- c(yg2,"H3K4me3")
      ys <- c(ys,y.pos)
      
      y.tmp3 <- y.pos
      m.temp3 <- m
      mark.tmp3 <- mark
      s <- min(c(min(c(mark.tmp3$start,mark.tmp3$end)),s)) - 100
      e <- max(c(max(c(mark.tmp3$start,mark.tmp3$end)),e)) + 100
      
      g2 <- g2 + geom_segment(data=mark.tmp3,aes(x=start,xend=end,y=y.tmp3,yend=y.tmp3,colour=cond),alpha=0.5,size=5) + geom_label(aes(x=s,y=y.tmp3,label=m.temp3))
      y.pos <- y.pos + 0.15
    }
    
    k4me3 <- reduce(GSE116695.peaks$K4me3_1,)
    k4me3$cond <- "5h"
    k4me3_naive <- reduce(GSE116695.peaks$Naive_H3K4me3)
    k4me3_naive$cond <- "0h"
    k4me3_24 <- reduce(GSE116695.peaks$E24h_K4me3)
    k4me3_24$cond <- "24h"
    k4me3 <- c(k4me3,k4me3_naive,k4me3_24)
    k4me3 <- data.table(as.data.frame(k4me3[as.vector(seqnames(k4me3))==chrom & ((start(k4me3)>s & start(k4me3)<e) | (end(k4me3)>s & end(k4me3)<e))]))
    
    if (nrow(k4me3)>0){
      yg2 <- c(yg2,"k4me3")
      ys <- c(ys,y.pos)
      
      y.tmp4me3 <- y.pos
      s <- min(c(min(c(k4me3$start,k4me3$end)),s)) - 100
      e <- max(c(max(c(k4me3$start,k4me3$end)),e)) + 100
      
      g2 <- g2 + geom_segment(data=k4me3,aes(x=start,xend=end,y=y.tmp4me3,yend=y.tmp4me3,color=cond),alpha=0.5,size=5) + geom_label(aes(x=s,y=y.tmp4me3,label="k4me3"))
      y.pos <- y.pos + 0.15
    }
    
    m <- "H3K27me3"; mark <-   format.mark(H3,m,s,e,chrom)
    if (nrow(mark)>0){
      yg2 <- c(yg2,"H3K27me3")
      ys <- c(ys,y.pos)
      
      y.tmp <- y.pos
      m.temp <- m
      mark.tmp <- mark
      s <- min(c(min(c(mark.tmp$start,mark.tmp$end)),s)) - 100
      e <- max(c(max(c(mark.tmp$start,mark.tmp$end)),e)) + 100
      
      g2 <- g2 + geom_segment(data=mark.tmp,aes(x=start,xend=end,y=y.tmp,yend=y.tmp,colour=cond),alpha=0.5,size=5) + geom_label(aes(x=s,y=y.tmp,label=m.temp))
      #y.pos <- y.pos + 0.1
    }
    
    polRA_Th1 <- reduce(c(GSE62486.peaks$Th1_reactivated_PolII,GSE62486.peaks$Th1_RA_PolII))
    polRA_Th1$cond <- "Th1"
    polRA_Th2 <- reduce(c(GSE62486.peaks$Th2_reactivated_PolII,GSE62486.peaks$Th2_RA_PolII))
    polRA_Th2$cond <- "Th2"
    RA_PolII <- c(polRA_Th1,polRA_Th2)
    RA_PolII <- data.table(as.data.frame(RA_PolII[as.vector(seqnames(RA_PolII))==chrom & ((start(RA_PolII)>s & start(RA_PolII)<e) | (end(RA_PolII)>s & end(RA_PolII)<e))]))
    
    if (nrow(RA_PolII)>0){
      #yg2 <- c(yg2,"RA_PolII")
      #ys <- c(ys,y.pos)
      
      y.tmpRApol <- 0.2
      s <- min(c(min(c(RA_PolII$start,RA_PolII$end)),s)) - 100
      e <- max(c(max(c(RA_PolII$start,RA_PolII$end)),e)) + 100
      
      g <- g + geom_segment(data=RA_PolII,aes(x=start,xend=end,y=y.tmpRApol,yend=y.tmpRApol,color=cond),alpha=0.5,size=5) + geom_label(aes(x=s,y=y.tmpRApol,label="RA_PolII"))
      #y.pos <- y.pos + 0.15
    }
    
    PolUS_Th1 <- reduce(c(GSE62486.peaks$Th1_unstimulated_PolII,GSE62486.peaks$Th1_US_PolII))
    PolUS_Th1$cond <- "Th1"
    PolUS_Th2 <- reduce(c(GSE62486.peaks$Th2_unstimulated_PolII,GSE62486.peaks$Th2_US_PolII))
    PolUS_Th2$cond <- "Th2"
    US_PolII <- c(PolUS_Th1,PolUS_Th2)
    US_PolII <- data.table(as.data.frame(US_PolII[as.vector(seqnames(US_PolII))==chrom & ((start(US_PolII)>s & start(US_PolII)<e) | (end(US_PolII)>s & end(US_PolII)<e))]))
    
    if (nrow(US_PolII)>0){
      #yg2 <- c(yg2,"US_PolII")
      #ys <- c(ys,y.pos)
      
      y.tmpUSpol <- 0.4
      s <- min(c(min(c(US_PolII$start,US_PolII$end)),s)) - 100
      e <- max(c(max(c(US_PolII$start,US_PolII$end)),e)) + 100
      
      g <- g + geom_segment(data=US_PolII,aes(x=start,xend=end,y=y.tmpUSpol,yend=y.tmpUSpol,color=cond),alpha=0.5,size=5) + geom_label(aes(x=s,y=y.tmpUSpol,label="US_PolII"))
      #y.pos <- y.pos + 0.15
    }
    
    g <- g + ylim(-1,1.3) + xlim(s,e)
    g2 <- g2 + ylim(1.2,y.pos+0.1) + xlim(s,e) + theme_bw() +
      scale_colour_manual(values=c(DMSO=colours()[124],BAY=colours()[184],RA=colours()[34],Th1=colours()[54],Th2=colours()[613],"0h"=colours()[124],"5h"=colours()[457],"24h"=colours()[144],up="red",down=colours()[131],stable="black",no.data="darkgrey",opening="red",closing=colours()[124],G0=colours()[124],G1="red")) 
    if (nrow(ATAC.select)>0){
      g2 <- g2 + geom_segment(data=ATAC.select,aes(x=center,xend=center,y=1.2,yend=y.pos,colour=status),size=0.1,linetype="dashed")
    }
    list(g=g,g2=g2,yg2=yg2,ys=ys,coords=paste0("The displayed region covers ",chrom,":",s,"-",e, "(",round((e-s)/1000000,3),"Mb)"),chrom=chrom,s=s,e=e)
  })
  
  
  output$chipSeq <- renderPlot({
    plot.TAD.summ()$g2
  })
  output$CHiC.TAD.plot <- renderPlot({
    plot.TAD.summ()$g
  })
  
  
  
  zoomin.CHiC <- reactive({
    g <-  plot.TAD.summ()$g
    coord <- ranges2$x
    if (is.null(coord[1])){
      coord <- c(1,2)
    }
    g <- g + coord_cartesian(xlim = coord,expand = F)
    g2 <-  plot.TAD.summ()$g2
    g2 <- g2 + scale_x_continuous(limits = coord,expand = c(0,0)) + geom_label(data=data.frame(x=rep(coord[1],length(plot.TAD.summ()$ys)),end=0,y=plot.TAD.summ()$ys,label=plot.TAD.summ()$yg2),aes(x=x,y=y,label=label))
    
    list(g=g,g2=g2,coords=paste0("The zoomed region covers ",plot.TAD.summ()$chrom,":",round(coord[1]),"-",round(coord[2]), "(",round((coord[2]-coord[1])/1000000,3),"Mb)"),chrom=plot.TAD.summ()$chrom,s=coord[1],e=coord[2])
    
  })
  
  output$zoomin.chipSeq <- renderPlot({
    zoomin.CHiC()$g2
  })
  
  output$zoomin.CHiC <- renderPlot({
    zoomin.CHiC()$g
  })
  
  output$CHiC.TAD.plot.coords <- renderText({
    plot.TAD.summ()$coords
  })
  
  output$CHiC.TAD.plot.coords.zoom <- renderText({
    zoomin.CHiC()$coords
  })
  
  hover.info <- reactive({
    homer.tmp$db.clust <- as.vector(homer.tmp$db.clust)
    homer.tmp$db.clust[!is.na(homer.tmp$cl1.split)] <- homer.tmp$cl1.split[!is.na(homer.tmp$cl1.split)] 
    homer.num$Interaction_Number <- NULL
    homer.num <- merge(external.chip.ATAC,homer.num,by="ATAC.id")
    colnames(homer.num) <- unlist(sapply(strsplit(colnames(homer.num),"[.][.]"),"[[",1))
    homer.num <- as.data.frame(homer.num)
    zoomin <- zoomin.CHiC()
    hover.coords <- input$plot_hover
    if (!is.null(hover.coords)){
      hover.coords.GR <- GRanges(zoomin$chrom,IRanges(hover.coords$x,hover.coords$x))
    }
    else {
      hover.coords$y <- 0
      hover.coords.GR <- GRanges(zoomin$chrom,IRanges(0,0))
    }
    
    res <- data.table(message="no element selected")
    if (hover.coords$y>1.15 & hover.coords$y<1.25){
      peaks.ATAC$center <- (start(peaks.ATAC)+end(peaks.ATAC))/2
      start(peaks.ATAC) <- start(peaks.ATAC) - 200
      end(peaks.ATAC) <- end(peaks.ATAC) + 200
      ovl <- findOverlaps(peaks.ATAC,hover.coords.GR,select="first")
      if (sum(!is.na(ovl))>0){
        ATAC <- data.table(as.data.frame(peaks.ATAC[!is.na(ovl)]))
        ATAC <- ATAC[which.min(abs(center-start(hover.coords.GR)))]
        id <- ATAC$id
        db.clust <- homer.tmp$db.clust[match(id,homer.tmp$ATAC.id)]
        TFs.num <- as.numeric(homer.num[match(id,homer.num$ATAC.id),-1])
        TFs <- colnames(homer.num)[-1]
        TFs <- TFs[TFs.num>0]
        TFs.num <- TFs.num[TFs.num>0]
        TFs <- paste(paste(TFs,TFs.num,sep=":"),collapse=", ")
        res <- ATAC[,list(id,tSNE=db.clust,status,CG,CpG,genomic,SYMBOL,RNAseq=expr,time.course=time.course$time.course[match(SYMBOL,time.course$SYMBOL)],H3K27me3,H3K4me1,H3K4me3,H3K27Ac,TFs=TFs)]
      }
    }
    else if (hover.coords$y>0.95 & hover.coords$y<1.05){
      frag.stats <- per.fragment.interactions.total()
      ovl <- findOverlaps(all.fragments,hover.coords.GR,select="first")
      if (sum(!is.na(ovl))>0){
        tmp <- all.fragments[!is.na(ovl)]
        id <- ifelse(tmp$is.bait=="other",tmp$IDs,tmp$gene)
        res <- frag.stats[name%in%id]
        res$pos <- res$start <- res$chr <- NULL
        res$in.genebody <- tmp$gene.bodies
        res$RNAseq <- ifelse(tmp$is.bait!="other",RNAseq[match(tmp$gene,ID)]$status,ifelse(res$in.genebody=="no.gene","intergenic",RNAseq[match(res$in.genebody,ID)]$status))
        res$time.course <- ifelse(tmp$is.bait!="other",time.course[match(tmp$gene,SYMBOL)]$time.course,ifelse(res$in.genebody=="no.gene","intergenic",time.course[match(res$in.genebody,SYMBOL)]$time.course))
      }
      
    }
    else if (hover.coords$y>-0.05 & hover.coords$y<0.05){
      ovl <- findOverlaps(refSeq,hover.coords.GR,select="first")
      if (sum(!is.na(ovl))>0){
        tmp <- data.table(as.data.frame(refSeq[!is.na(ovl)]))
        res <- tmp[,list(SYMBOL,strand,RNAseq=RNAseq[match(SYMBOL,ID)]$status,time.course=time.course$time.course[match(SYMBOL,time.course$SYMBOL)])]
        
      }
    }
    return(res)
  })
  
  output$hover.info <- renderDataTable({
    hover.info()
  })
  
  interaction.network.CHiC <- eventReactive(input$go.net.CHiC,{
    SYMBOL <- input$SYMBOL.net.CHiC
    net.type <- input$type.net.CHiC
    diff.interactions <- all.interactions[extend.classif!="stable"]
    diff.interactions$idx.int <- 1:nrow(diff.interactions)
    selected.interactions <- diff.interactions[name1==SYMBOL | name2==SYMBOL]
    switch(net.type,
           only={
             selected.interactions 
           },
           extended={
             all.ids <- unique(c(selected.interactions$name1,selected.interactions$name2))
             temp <- diff.interactions[name1%in%all.ids | name2%in%all.ids]
             inter.sel <- selected.interactions$idx.int
             while (sum(!temp$idx.int%in%inter.sel)>0){
               inter.sel <- unique(c(inter.sel,temp$idx.int))
               all.ids <- unique(c(temp$name1,temp$name2))
               temp <- diff.interactions[name1%in%all.ids | name2%in%all.ids]
               
             }
             selected.interactions <- diff.interactions[name1%in%all.ids & name2%in%all.ids]
           }
           )
    selected.interactions
  })
  
  output$interaction.network.CHiC <- renderPlotly({
    g <- make.network.interaction(interaction.network.CHiC(),RNAseq)
  })
  
  output$interaction.network.CHiC.tab <- renderDataTable({
    interaction.network.CHiC()[,list(chr1,locus1,locus2,name1,name2,type,differential=extend.classif,distance)]
  })
  
  output$CHiC.TAD.plot.coords.network <- renderText({
    tab <- interaction.network.CHiC()
    chrom <- unique(tab$chr1)
    s <- min(c(tab$locus1,tab$locus2))
    e <- max(c(tab$locus1,tab$locus2))
    paste0("The displayed region covers ",chrom,":",s,"-",e, "(",round((e-s)/1000000,3),"Mb)")
  })
  
  
  

  #################################################################################################
  ## integration - Figure 2c Burren
  #################################################################################################
  
  ATAC.bait.RNA <- reactive({
    message("Preparing plot from Burren et al. ...")
    
    # all.fragments
    # peaks.ATAC
    # all.interactions
    # RNAseq
    ovl <- findOverlaps(peaks.ATAC,all.fragments[all.fragments$is.bait=="other"])
    ATAC.other <- cbind(as.data.frame(mcols(peaks.ATAC))[queryHits(ovl),c("id","logFC","status","genomic")],as.data.frame(mcols(all.fragments[all.fragments$is.bait=="other"]))[subjectHits(ovl),c("IDs")])
    ATAC.other <- data.table(ATAC.other)
    colnames(ATAC.other) <- c("ATAC.id","ATAC.logFC","ATAC.status","ATAC.genomic","frag.id")
    tmp.interactions <- all.interactions[type=="bo",list(bait=name1,other=name2,CHiC.status=extend.classif,CHiC.distance=distance)]
    tmp.interactions$RNAseq.status <- RNAseq[match(tmp.interactions$bait,RNAseq$ID)]$status
    tmp.interactions$RNAseq.logFC <- RNAseq[match(tmp.interactions$bait,RNAseq$ID)]$logFC
    m <- matches(ATAC.other$frag.id,tmp.interactions$other,all.x=F)
    
    ATAC.other.CHiC <- cbind(ATAC.other[m$x,],tmp.interactions[m$y,])
    ATAC.other.CHiC <- ATAC.other.CHiC[!is.na(RNAseq.logFC)]
    ATAC.other.CHiC[is.na(ATAC.other.CHiC)] <- "no"
    ATAC.other.CHiC <- ATAC.other.CHiC[,list(nb.frag=length(CHiC.distance)),by=c("ATAC.id","ATAC.status","CHiC.status","bait","RNAseq.logFC","RNAseq.status","other")]
    ATAC.other.CHiC <- ATAC.other.CHiC[,list(nb.frag=length(other)),by=c("ATAC.status","CHiC.status","bait","RNAseq.logFC","RNAseq.status")]
    ATAC.other.CHiC$group <- paste(ATAC.other.CHiC$CHiC.status,ATAC.other.CHiC$ATAC.status,sep="_")
    
    ATAC.other.CHiC.sn <- ATAC.other.CHiC[!duplicated(bait)][grep("stable_",group)]
    ATAC.other.CHiC <- rbind(ATAC.other.CHiC[!grepl("stable_",group)],ATAC.other.CHiC.sn)
    
    g1 <- ggplot(ATAC.other.CHiC,aes(x=group)) + geom_bar(fill="lightgrey",color="darkgrey") + scale_y_log10() + theme_bw()
    ATAC.other.CHiC.stats <- ATAC.other.CHiC[,list(average=mean(RNAseq.logFC),CI=qnorm(0.975)*sd(RNAseq.logFC)/sqrt(length(RNAseq.status))),by=c("group","ATAC.status","CHiC.status")]
    g2 <- ggplot(ATAC.other.CHiC.stats,aes(x=group,y=average)) + ylab("bait logFC expression +/- 95%CI") +
      geom_hline(yintercept = 0,linetype="dashed",size=0.3) + geom_point(aes(color=CHiC.status))  + theme_bw() + theme(legend.position="none") +
      geom_errorbar(aes(ymin=average-CI, ymax=average+CI,color=CHiC.status), width=0,size=0.5) + scale_color_manual(values=c(gain="red",loss=colours()[124],stable="darkgrey"))
   
    multiplot(g1,g2,cols=1)
      })
  
  output$ATAC.bait.RNA <- renderPlot({
    ATAC.bait.RNA()
  })
  
  #################################################################################################
  ## integration - tSNE
  #################################################################################################
  
  
  prepare.shaps <- reactive({
    withProgress(message="Preparing data for SHAP analysis",{
      diff.fragments <- all.interactions[extend.classif!="stable",list(id1,id2,name1,name2,type1,type2,type,diff.int=extend.classif)]
      
      m <- match(diff.fragments$name1,RNAseq$ID)
      diff.fragments$expr1 <- RNAseq$status[m]
      diff.fragments$logFC1 <- RNAseq$logFC[m]
      
      m <- match(diff.fragments$name2,RNAseq$ID)
      diff.fragments$expr2 <- RNAseq$status[m]
      diff.fragments$logFC2 <- RNAseq$logFC[m]
      
      diff.fragments <- diff.fragments[!is.na(expr1)]
      diff.fragments <- diff.fragments[type=="bo" | (type=="bb" & !is.na(expr2))]
      diff.fragments <- rbind(diff.fragments[,list(id=id1,name=name1,type=type1,B2B=type,diff.int,diff.expr=expr1,logFC.expr=logFC1)],diff.fragments[,list(id=id2,name=name2,type=type2,B2B=type,diff.int,diff.expr=ifelse(type=="bo",expr1,expr2),logFC.expr=ifelse(type=="bo",logFC1,logFC2))])
      diff.fragments[,logFC.expr:=NULL]

      diff.fragments <- unique(diff.fragments[,-4,with=F],by=NULL)
      
      diff.fragments$seqnames <- as.vector(seqnames(all.fragments))[match(diff.fragments$id,all.fragments$IDs)]
      diff.fragments$start <- as.vector(start(all.fragments))[match(diff.fragments$id,all.fragments$IDs)]
      diff.fragments$end <- as.vector(end(all.fragments))[match(diff.fragments$id,all.fragments$IDs)]
      
      diff.fragments <- as(diff.fragments[,list(seqnames,start,end,id=id,diff.int,bait.SYMBOL=name,type,bait.status=diff.expr)],"GRanges")
      ovl <- findOverlaps(peaks.ATAC,diff.fragments)
      nb.int.diff <- length(unique(subjectHits(ovl)))
      peaks.ATAC.frag <- peaks.ATAC[queryHits(ovl)]
      peaks.ATAC.frag$in.diff.frag <- "yes" 
      peaks.ATAC.frag$diff.int <- diff.fragments$diff.int[subjectHits(ovl)]
      peaks.ATAC.frag$id.frag <- diff.fragments$id[subjectHits(ovl)]
      peaks.ATAC.frag$bait.symbol<- diff.fragments$bait.SYMBOL[subjectHits(ovl)]
      peaks.ATAC.frag$bait.status<- diff.fragments$bait.status[subjectHits(ovl)]
      peaks.ATAC.frag$frag.type<- diff.fragments$type[subjectHits(ovl)]
      peaks.ATAC.frag$all.marks <- paste(peaks.ATAC.frag$H3K27me3,"/",peaks.ATAC.frag$H3K4me3,peaks.ATAC.frag$H3K4me1,peaks.ATAC.frag$H3K27Ac,sep="/") 
      m <- match(peaks.ATAC.frag$id,homer.num$ATAC.id)
      homer.num <- homer.num[m,]
      homer.num$Feature <- paste(peaks.ATAC.frag$diff.int,peaks.ATAC.frag$bait.status,sep="_")
      homer.num$Interaction_Number <- 1:nrow(homer.num)
      #homer.num <- homer.num[order(homer.num$Interaction_Number),]
      homer.num$ATAC.status <- peaks.ATAC.frag$status[!is.na(m)]
      homer.num$frag.type <- peaks.ATAC.frag$frag.type[!is.na(m)]
      homer.num$genomic <- peaks.ATAC.frag$genomic[!is.na(m)]
      homer.num$CG <- peaks.ATAC.frag$CG[!is.na(m)]
      homer.num$bait.status <- peaks.ATAC.frag$bait.status[!is.na(m)]
      homer.num$H3K4me1 <- peaks.ATAC.frag$H3K4me1[!is.na(m)]
      homer.num$H3K4me3 <- peaks.ATAC.frag$H3K4me3[!is.na(m)]
      homer.num$H3K27Ac <- peaks.ATAC.frag$H3K27Ac[!is.na(m)]
      homer.num$H3K27me3 <- peaks.ATAC.frag$H3K27me3[!is.na(m)]
      
      homer.num$all.marks <- peaks.ATAC.frag$all.marks[!is.na(m)]
      
      homer.num$SYMBOL.ATAC <- peaks.ATAC.frag$SYMBOL[!is.na(m)]
      homer.num$frag.ID <- peaks.ATAC.frag$bait.symbol[!is.na(m)]
      
      homer.num <- cbind(homer.num,external.chip.ATAC[match(homer.num$ATAC.id,external.chip.ATAC$ATAC.id),!colnames(external.chip.ATAC)%in%c("ATAC.id","db.clust"),with=F])
      
      
      shaps <- na.omit(shaps)
      shaps <- shaps[!grepl("stable",shaps$feature),]
      
      m <- match(shaps$ID,homer.num$Interaction_Number)
      
      homer.tmp <- homer.num[na.omit(m),!duplicated(colnames(homer.num))]
      row.names(homer.tmp) <- homer.tmp$Interaction_Number
      homer.tmp$db.clust <- factor(tsne.group$db.clust)
      homer.tmp$cl1.split <- tsne.group$cl1.split
      
      colnames(homer.tmp) <- unlist(sapply(strsplit(colnames(homer.tmp),"[.][.]"),"[[",1))
      to.discard <- which(apply(homer.tmp[,1:436],2,sum)==0)
      homer.tmp <- homer.tmp[,-to.discard]
      return(list(homer.tmp=homer.tmp,shaps=shaps,tsne.group=tsne.group,tsne.center=tsne.center,peaks.ATAC.frag=peaks.ATAC.frag,diff.fragments=diff.fragments))
      
    })
  })
  
  output$TFs.tsne.list <- renderUI({
    homer.tmp <- prepare.shaps()$homer.tmp
    toselect <- colnames(homer.tmp)[!colnames(homer.tmp)%in%c("Interaction_Number","frag.ID","SYMBOL.ATAC","ATAC.id")]
    selectInput("select.TFs.tsne.list","Choose feature",choices=rev(toselect),selected="db.clust")
  })
  
  tSNE.shaps.icis <- eventReactive(input$go.tSNE.shaps ,{
    shaps.data <- prepare.shaps()
    homer.tmp <- shaps.data$homer.tmp
    shaps <- shaps.data$shaps
    tsne.group <- shaps.data$tsne.group
    tsne.center <- shaps.data$tsne.center
    
    tsne.center$to.eval <- tsne.center$db.clust
    tsne.center$label <- tsne.center$db.clust
    
    symbols <- rbindlist(lapply(colnames(homer.tmp[1:416]),function(x){
      x <- sub("T[.]box","Tbox",x)
      tmp <- unlist(strsplit(x,"[.]"))
      TF <- paste(tmp[-length(tmp)],collapse=".")
      type <- tmp[length(tmp)]
      data.table(tf=TF,type)
    }))
    TF.labels <- apply(homer.tmp[,1:416],1,function(x,symbols){
      paste(symbols[x>0],collapse ="\n")
    },symbols$tf)
    tsne.group$label <- paste(homer.tmp$Feature,paste0("genomic:",homer.tmp$genomic,":",homer.tmp$SYMBOL.ATAC),paste0("fragment:",homer.tmp$frag.ID),paste0("H3K4me3:",homer.tmp$H3K4me3),paste0("H3K4me1:",homer.tmp$H3K4me1),paste0("H3K27Ac:",homer.tmp$H3K27Ac),TF.labels,sep="\n")
    tsne.group$to.eval <- homer.tmp[,input$select.TFs.tsne.list]
    color.clust <- c(rep(c(brewer.pal(8, "Set1"),brewer.pal(8, "Set2"),brewer.pal(8, "Accent"),brewer.pal(8, "Spectral"),brewer.pal(8, "Pastel1"),brewer.pal(10, "BrBG")),3))
    names(color.clust) <- as.character(1:length(color.clust))
    if (is.numeric(tsne.group$to.eval)){
      tsne.group <- tsne.group[order(tsne.group$to.eval,decreasing=F),]
    }
    g <- ggplot(aes(x=tSNE1, y=tSNE2, group=to.eval,label=label), data=tsne.group) + 
      geom_point(aes(x=tSNE1, y=tSNE2, colour=to.eval),size=0.5) +   
      theme(panel.background = element_rect(fill = "white", colour = "black")) +
      geom_hline(yintercept=0,linetype="dashed",size=0.3) +
      geom_vline(xintercept=0,linetype="dashed",size=0.3) 
    if (is.numeric(tsne.group$to.eval)){
      g <- g + scale_colour_gradientn(colours = colours()[c(234,420,34,36,24)])
        #scale_colour_gradient2(low = "blue",mid = "lightgrey",high = "red",midpoint = 0)
      #colours = colours()[c(234,420,34,36,24)]
    } 
    else if (input$select.TFs.tsne.list%in%c("H3K4me1","H3K4me3","H3K27Ac","H3K27me3")) {
      g <- g + scale_colour_manual(values=c("00"="grey","10"="blue","01"="red","11"="green"))
    }
    else if (input$select.TFs.tsne.list=="all.marks") {
      names(color.clust) <- unique(tsne.group$to.eval)
      g <- g + scale_colour_manual(values=color.clust)
    }
    else if (input$select.TFs.tsne.list=="cl1.split") {
      names(color.clust) <- unique(na.omit(unique(tsne.group$cl1.split)))
      g <- g + scale_colour_manual(values=color.clust)
    }
    else {
      g <- g + scale_colour_manual(values=c(intergenic="grey",tss=colours()[614],prom=colours()[76],gene.body=colours()[613],closing=colours()[131],opening="red","stable"="grey",gain_up="red",gain_stable=colours()[235],gain_down=colours()[53],loss_up=colours()[96],loss_stable=colours()[613],loss_down=colours()[124],up="red",down=colours()[124],gain="red",loss=colours()[124],b=colours()[613],o=colours()[34],color.clust)) 
    }
    if (input$show.db.clust){
      g <- g + geom_text(data=tsne.center)
    }
    g
  })
  
  output$tSNE.shaps.icis <- renderPlot({
    tSNE.shaps.icis()
  })
  
  
  brush_info <- reactive({
    data.brush <- input$plot_brush
    tsne.group$db.clust <- as.vector(tsne.group$db.clust)
    tsne.group$db.clust[!is.na(tsne.group$cl1.split)] <- tsne.group$cl1.split[!is.na(tsne.group$cl1.split)] 
    homer.tmp <- prepare.shaps()$homer.tmp
    homer.tmp$db.clust <- as.vector(homer.tmp$db.clust)
    homer.tmp$db.clust[!is.na(homer.tmp$cl1.split)] <- homer.tmp$cl1.split[!is.na(homer.tmp$cl1.split)] 
    
    
    select.brush <- tsne.group$tSNE1>=data.brush$xmin & tsne.group$tSNE1<=data.brush$xmax & tsne.group$tSNE2>=data.brush$ymin & tsne.group$tSNE2<=data.brush$ymax
    enriched <- homer.tmp[select.brush,c(1:416,432:466)]
    annots <- homer.tmp[select.brush,c("ATAC.id","db.clust","Feature","ATAC.status","frag.type","genomic","SYMBOL.ATAC","H3K4me1","H3K4me3","H3K27Ac","H3K27me3")]
    #annots$nb.annots <- apply(enriched,1,function(x) sum(x>0))
    annots.rows <- annots
    tsne.group <- tsne.group[select.brush,]
    if (input$show.all.points.in.area){
      nb.f <- table(tsne.group$db.clust)
      cl <- names(nb.f)[which.max(nb.f)]
      enriched <- enriched[annots.rows$db.clust==cl,]
      annots.rows <- annots.rows[annots.rows$db.clust==cl,]
    }
    
    list(enriched=enriched,annots.rows=annots.rows)
  })
  
  select.features <- reactive({
    data <- brush_info()
    enriched <- data$enriched
    annots.rows <- data$annots.rows
    
    if (input$show.feature.in.area){
      enriched <- cbind(enriched,annots.rows[,!colnames(annots.rows)%in%colnames(enriched),])
      
      if(is.numeric(enriched[,input$select.TFs.tsne.list])){
        sliderInput("nb.min.feature","Choose minimal number of feature per peak to select",min = 0,max=max(enriched[,input$select.TFs.tsne.list]),value = 0,step=1)
      }
      else if(!is.numeric(enriched[,input$select.TFs.tsne.list])){
        choices <- unique(as.vector(enriched[,input$select.TFs.tsne.list]))
        checkboxGroupInput("select.features.values","Select feature(s) to keep peaks carrying it/them.",choices = choices,inline = T)
      }
    }
  })
  
  output$select.features <- renderUI({
    select.features()
  })
  
  brush_info2 <- reactive({
    data <- brush_info()
    enriched <- data$enriched
    annots.rows <- data$annots.rows
    
    if (input$show.feature.in.area){
      n <- ncol(enriched)
      enriched <- cbind(enriched,annots.rows[,!colnames(annots.rows)%in%colnames(enriched),])
      if(is.numeric(enriched[,input$select.TFs.tsne.list])){
        annots.rows <- annots.rows[enriched[,input$select.TFs.tsne.list]>=input$nb.min.feature,]
        enriched <- enriched[enriched[,input$select.TFs.tsne.list]>=input$nb.min.feature,]
      }
      if(!is.numeric(enriched[,input$select.TFs.tsne.list])){
        annots.rows <- annots.rows[enriched[,input$select.TFs.tsne.list]%in%input$select.features.values,]
        enriched <- enriched[enriched[,input$select.TFs.tsne.list]%in%input$select.features.values,]
      }
      enriched <- enriched[,1:n]
    }
    list(enriched=enriched,annots.rows=annots.rows)
  })
  
  output$nb.frag.selected <- renderText({
    text <- brush_info2()$annots.rows
    text <- paste(nrow(text),"ATAC peaks are selected.")
  })
  
  heatmap.shaps.clust <- reactive({
    brush.data <- brush_info2()
    enriched <- brush.data$enriched
    annots.rows <- brush.data$annots.rows
    annots.rows <- annots.rows[,c("ATAC.id","db.clust","Feature","ATAC.status","SYMBOL.ATAC","frag.type","genomic","H3K4me1","H3K4me3","H3K27Ac","H3K27me3")]
    enriched <- enriched[!duplicated(annots.rows$ATAC.id),]
    annots.rows <- annots.rows[!duplicated(annots.rows$ATAC.id),]
    row.names(enriched) <- annots.rows$ATAC.id
    row.names(annots.rows) <- annots.rows$ATAC.id
    annots.rows <- annots.rows[,-1]
    enriched <- enriched[,apply(enriched,2,function(x) (sum(x>0)/nrow(enriched))>input$min.pc.feature)]
    
    annots.rows <- annots.rows[row.names(annots.rows)%in%row.names(enriched),]
    enriched <- enriched[,order(apply(enriched,2,function(x) sum(x>0)/nrow(enriched)),decreasing=T)]
    #annots.rows <- annots.rows[,-1]
    clust <-  pheatmap(t(enriched),
                       clustering_distance_cols = "binary",
                       clustering_method = "ward.D",silent=T)
    nb.clust <- input$nb.clust.heatmap
    cluster=factor(cutree(clust$tree_col, k = nb.clust))
    annots.rows$heatmap <- cluster[match(row.names(annots.rows),names(cluster))]
    annots.rows <- annots.rows[match(names(cluster),row.names(annots.rows)),]
    #annots.rows$position <- 1:nrow(annots.rows)
    
    data <- list(enriched=enriched,annots.rows=annots.rows,clust.heatmap=clust)    
  })
  
  heatmap.shaps.clust.plot <- reactive({
    data <- heatmap.shaps.clust()
    data$annots.rows$SYMBOL.ATAC <- NULL
    homer.tmp <- prepare.shaps()$homer.tmp
    homer.tmp$db.clust <- as.vector(homer.tmp$db.clust)
    homer.tmp$db.clust[!is.na(homer.tmp$cl1.split)] <- homer.tmp$cl1.split[!is.na(homer.tmp$cl1.split)] 
    color.clust <- c(rep(c(brewer.pal(8, "Set1"),brewer.pal(8, "Set2"),brewer.pal(8, "Accent"),brewer.pal(8, "Spectral"),brewer.pal(8, "Pastel1"),brewer.pal(10, "BrBG")),3))
    names(color.clust) <- c(unique(homer.tmp$db.clust))
    color.clust <- color.clust[!is.na(names(color.clust))]
    pheatmap(t(data$enriched),
             clustering_distance_cols = "binary",
             #clustering_distance_rows = "binary",
             cluster_rows = F,
             color=colorRampPalette(colours()[c(1,34,35,24)])(100),
             show_rownames = T,show_colnames = F,
             border_color = NA, annotation_col = data$annots.rows,
             clustering_method = "ward.D",
             annotation_colors = list(db.clust=color.clust[names(color.clust)%in%data$annots.rows$db.clust],genomic=c(intergenic="grey",tss=colours()[614],prom=colours()[76],gene.body=colours()[613]),frag.type=c(o="black",b="grey"),Feature=c(gain_up="red",gain_down=colours()[53],loss_up=colours()[96],loss_down=colours()[124]),ATAC.status=c(closing=colours()[131],opening="red",stable="grey"),H3K4me1=c("00"="grey","10"="blue","01"="red","11"="green"),H3K4me3=c("00"="grey","10"="blue","01"="red","11"="green"),H3K27me3=c("00"="grey","10"="blue","01"="red","11"="green"),H3K27Ac=c("00"="grey","10"="blue","01"="red","11"="green")))
  })
  
  output$heatmap.shaps.clust.plot <- renderPlot({
    heatmap.shaps.clust.plot()
  })
  
  output$feature.shaps.clust <- renderPlot({
    data <- heatmap.shaps.clust()
    enriched <- data$enriched
    features <- sapply(enriched,function(x) sum(x>0)/nrow(enriched))
    features <- data.frame(features=names(features),freq=as.numeric(features))
    features$features <- factor(as.vector(features$features),levels=features$features[order(features$freq,decreasing = F)])
    ggplot(features,aes(x=features,y=freq)) + geom_col(fill="grey",colour="black") + theme_bw() +
      #theme(axis.text.x = element_text(angle = 90, hjust = 1,size=12)) + 
      coord_flip()
  })
  
  
  output$brush_info <- renderDataTable({
    peaks.ATAC.frag <- prepare.shaps()$peaks.ATAC.frag
    annots.rows <- heatmap.shaps.clust()$annots.rows
    annots.rows$RNAseq.ATAC <- RNAseq[match(annots.rows$SYMBOL.ATAC,ID)]$status
    annots.rows$cinetic.ATAC <- time.course[match(annots.rows$SYMBOL.ATAC,SYMBOL)]$time.course
    m <- match(row.names(annots.rows),peaks.ATAC.frag$id)
    annots.rows$frag.id <-  peaks.ATAC.frag$bait.symbol[m]
    annots.rows$RNAseq.bait <- RNAseq[match(annots.rows$frag.id,ID)]$status
    annots.rows$cinetic.bait <- time.course[match(annots.rows$frag.id,SYMBOL)]$time.course
    tSNE.info <- interacting.tSNE()
    m <- data.table(as.data.frame(matches(row.names(annots.rows),tSNE.info$ATAC.id1,all.x=T,all.y=F)))
    m$frag.id2 <- tSNE.info$fragment.id2[m$y]
    m <- m[,list(frag.id2=paste(unique(frag.id2),collapse=", ")),by=x]
    annots.rows$int.frag[m$x] <- m$frag.id2
    annots.rows[,!grepl("H3K",colnames(annots.rows))]
  })
  
  test.tab <- eventReactive(input$go.save.shaps.tsne,{
    withProgress(message = 'Dowloading files ...', {
      annots.rows <- heatmap.shaps.clust()$annots.rows
      annots.rows$RNAseq <- RNAseq[match(annots.rows$SYMBOL.ATAC,ID)]$status
      annots.rows$time.course <- time.course[match(annots.rows$SYMBOL.ATAC,SYMBOL)]$time.course
      
      t2 <- heatmap.shaps.clust()$enriched
      tab <- cbind(annots.rows,t2)
      write.table(tab,file.path(data.dir,"tables",paste0(input$prefix.save.shaps.tsne,".tab")),sep="\t",quote=F)
    })
    text <- paste("File",paste0(input$prefix.save.shaps.tsne,".tab"),"has been written in the table folder.")
  })
  
  output$download.done <- renderText({
    test.tab()
  })
  
  output$tab.shaps.clust.plot <- renderPlot({
    annots.rows <- heatmap.shaps.clust()$annots.rows
    g1 <- ggplot(annots.rows,aes(x=Feature)) + geom_bar(aes(fill=Feature),color="black",alpha=0.5) + theme_bw() + facet_wrap(~heatmap) + scale_fill_manual(values=c(gain_up="red",gain_down=colours()[53],loss_up=colours()[96],loss_down=colours()[124])) + theme(legend.position = "none")# + theme(axis.text.x = element_text(angle = 90, hjust = 1))
    g2 <- ggplot(annots.rows,aes(x=genomic)) + geom_bar(aes(fill=genomic),color="black",alpha=0.5) + theme_bw() + facet_wrap(~heatmap) + scale_fill_manual(values=c(intergenic="grey",tss=colours()[614],prom=colours()[77],gene.body=colours()[613])) + theme(legend.position = "none")# + theme(axis.text.x = element_text(angle = 90, hjust = 1))
    g3 <- ggplot(annots.rows,aes(x=ATAC.status)) + geom_bar(aes(fill=ATAC.status),color="black",alpha=0.5) + theme_bw() + facet_wrap(~heatmap) + scale_fill_manual(values=c(stable="grey",opening="red",closing=colours()[124])) + theme(legend.position = "none")# + theme(axis.text.x = element_text(angle = 90, hjust = 1))
    g4 <- ggplot(annots.rows,aes(x=frag.type)) + geom_bar(aes(fill=frag.type),color="black",alpha=0.5) + theme_bw() + facet_wrap(~heatmap) + scale_fill_manual(values=c(b="grey",o="black")) + theme(legend.position = "none")# + theme(axis.text.x = element_text(angle = 90, hjust = 1))
    multiplot(g1,g2,g3,g4,cols=2)
  })
  
  output$tab.shaps.clust.plot2 <- renderPlot({
    annots.rows <- heatmap.shaps.clust()$annots.rows
    annots.rows <- data.table(annots.rows)[,list(Feature,ATAC.status, heatmap,frag.type,genomic, histones=paste(H3K27me3,H3K4me3,H3K4me1,H3K27Ac,sep="/" ))]
    ggplot(annots.rows,aes(x=histones)) + geom_bar(fill="grey",color="black") + theme_bw() + facet_wrap(~heatmap,ncol=1,scales="free_y") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle("H3K27me3/H3K4me3/H3K4me1/H3K27Ac")
  })
  
  output$tab.shaps.clust.plot3 <- renderPlot({
    annots.rows <- heatmap.shaps.clust()$annots.rows
    annots.rows <- data.table(annots.rows)[,list(heatmap,H3K27me3,H3K4me3,H3K4me1,H3K27Ac)]
    annots.rows <- melt(annots.rows,id.vars="heatmap")
    annots.rows$value <- factor(as.vector(annots.rows$value),levels=rev(c("00","10","01","11")))
    ggplot(annots.rows,aes(x=variable)) + geom_bar(aes(fill=value),color="black",alpha=0.5) + theme_bw() + scale_fill_manual(values=c("00"="white","10"=colours()[142],"01"=colours()[131],"11"=colours()[613])) + facet_wrap(~heatmap,scales="free_y") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle("H3K27me3/H3K4me3/H3K4me1/H3K27Ac")
  })
  
  interacting.tSNE <- reactive({
    heatmap.info <- heatmap.shaps.clust()$annots.rows
    brush.data <- brush_info2()
    enriched <- brush.data$enriched
    annots.rows <- brush.data$annots.rows
    homer.tmp <- prepare.shaps()$homer.tmp
    
    homer.tmp$db.clust <- as.vector(homer.tmp$db.clust)
    homer.tmp$db.clust[!is.na(homer.tmp$cl1.split)] <- homer.tmp$cl1.split[!is.na(homer.tmp$cl1.split)] 
    
    shaps.data <- prepare.shaps()
    #all.interactions
    peaks.ATAC.frag <- shaps.data$peaks.ATAC.frag
    diff.fragments <- shaps.data$diff.fragments
    
    tmp <- all.interactions[extend.classif!="stable"]
    
    tmp <- rbind(tmp[,list(id1,id2,type1,type2,name1,name2,type,CHiC=extend.classif)],
                 tmp[,list(id1=id2,id2=id1,type1=type2,type2=type1,name1=name2,name2=name1,type,CHiC=extend.classif)])
    
    m <- as.data.frame(matches(annots.rows$ATAC.id,peaks.ATAC.frag$id,all.y=F))
    annots.rows <- annots.rows[m$x,]
    annots.rows$id.frag <- peaks.ATAC.frag$id.frag[m$y]
    annots.rows$bait.symbol <- peaks.ATAC.frag$bait.symbol[m$y]
    annots.rows$bait.status <- peaks.ATAC.frag$bait.status[m$y]
    
    m1 <- as.data.frame(matches(annots.rows$id.frag,tmp$id1,all.y=F))
    m1$cl <- annots.rows$db.clust[m1$x]
    
    m1$ATAC1 <- annots.rows$ATAC.id[m1$x]
    m1$SYMBOL <- tmp$name1[m1$y]
    m1$Feature1 <- annots.rows$Feature[m1$x]
    
    m1$CHiC.diff <- tmp$CHiC[m1$y]
    m1$id <- tmp$id2[m1$y]
    m1$name2 <- tmp$name2[m1$y]
    m1$bait.status <- diff.fragments$bait.status[match(annots.rows$id.frag[m1$x],diff.fragments$id)]
    
    m1prime <- as.data.frame(matches(m1$id,peaks.ATAC.frag$id.frag,all.y=F))
    m1 <- m1[m1prime$x,]
    m1$ATAC <- peaks.ATAC.frag$id[m1prime$y]
    m1prime <- as.data.frame(matches(m1$id,diff.fragments$id,all.y=F))
    m1 <- m1[m1prime$x,]
    m1$Feature2 <- paste(diff.fragments$diff.int[m1prime$y],diff.fragments$bait.status[m1prime$y],sep="_")
    
    m1$clust <- homer.tmp$db.clust[match(m1$ATAC,homer.tmp$ATAC.id)]
    m1$clust[is.na(m1$ATAC)] <- "no.ATAC"
    on.bait.num <- (1:nrow(annots.rows))[-unique(m1[m1$Feature1==m1$Feature2,]$x)]
    m <- data.table(m1[m1$Feature1==m1$Feature2,])
    m <- m[,list(ATAC.id=ATAC1,cluster1=cl,SYMBOL,Feature=Feature1,interacting.frag=name2,ATAC,cluster2=clust)]
    on.bait <- annots.rows[on.bait.num,]
    on.bait$CHiC <- unlist(sapply(strsplit(on.bait$Feature,"_"),"[[",1))
    on.bait$CHiC <- ifelse(on.bait$CHiC=="gain","nb.gained","nb.lost")
    int.fragments.stats <- per.fragment.interactions.total()
    int.fragments.stats <- int.fragments.stats[match(on.bait$bait.symbol,name)]
    int.fragments.stats <- melt(int.fragments.stats,id.vars="name")
    on.bait$interacting.frag <- paste(on.bait$CHiC,int.fragments.stats$value[int.fragments.stats$variable==on.bait$CHiC],sep=": ")
    m <- rbind(m,data.table(on.bait)[,list(ATAC.id,cluster1=db.clust,SYMBOL=bait.symbol,Feature=bait.status,interacting.frag,ATAC="on.bait",cluster2=db.clust)])
    
    m[is.na(cluster2)]$cluster2 <- "not.in.tSNE"   
    #m[cluster2=="0"]$cluster2 <- "not.in.cluster"
    m$heatmap.info <- heatmap.info$heatmap[match(m$ATAC.id,row.names(heatmap.info))]
    colnames(m) <- c("ATAC.id1","cluster1","fragment.id1","class","fragment.id2","ATAC.id2","cluster2","heatmap.info")
    m <- m[,c("heatmap.info","class","cluster1","ATAC.id1","fragment.id1","ATAC.id2","fragment.id2","cluster2")]
    m <- m[!is.na(m$heatmap.info),]
    m <- m[m$class!="stable",]
    m <- data.table(unique(as.data.frame(m)))
  })
  
 
  
  output$tab.interacting.tSNE <- renderDataTable({
    interacting.tSNE()
  })
  
  
  
  
  
  output$all.shaps.clust <- renderPlotly({
    homer.tmp <- prepare.shaps()$homer.tmp
    homer.tmp$db.clust <- as.vector(homer.tmp$db.clust)
    homer.tmp$db.clust[!is.na(homer.tmp$cl1.split)] <- homer.tmp$cl1.split[!is.na(homer.tmp$cl1.split)] 
    alluvial <- interacting.tSNE()
    noATAC <- nrow(alluvial[cluster2=="no.ATAC"])
    alluvial <- alluvial[ATAC.id2!="on.bait" & cluster2!="no.ATAC"]
    color.clust <- c(rep(c(brewer.pal(8, "Set1"),brewer.pal(8, "Set2"),brewer.pal(8, "Accent"),brewer.pal(8, "Spectral"),brewer.pal(8, "Pastel1"),brewer.pal(10, "BrBG")),3))
    names(color.clust) <- c(unique(homer.tmp$db.clust),"not.in.tSNE")
    color.clust <- color.clust[!is.na(names(color.clust))]
    ggplot(alluvial,aes(x=cluster2)) + ggtitle(paste(noATAC, "interacting fragments do not carry ATAC peaks")) + geom_bar(aes(fill=cluster2),color="black") + xlab("cluster of interacting fragments") + scale_fill_manual(values=c("0"="white",color.clust)) + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  })
  
  interacting.heatmap <- reactive({
    data.heatmap <- heatmap.shaps.clust()
    data.interaction <- interacting.tSNE()
    data.interaction <- data.interaction[data.interaction$cluster2!="no.ATAC" & data.interaction$ATAC.id2!="on.bait" ,]
    id1 <- data.interaction$ATAC.id1
    id2 <- data.interaction$ATAC.id2
    #colnames(homer.num) <- unlist(sapply(strsplit(colnames(homer.num),"[.][.]"),"[[",1))
    ids <- paste(id1,id2,sep=":")

    homer.num$Interaction_Number <- NULL
    data.chip <- merge(external.chip.ATAC,homer.num,by="ATAC.id")

    tab1 <- as.data.frame(data.chip[data.chip$ATAC.id%in%id1,])
    row.names(tab1) <- tab1$ATAC.id
    tab1$ATAC.id <- NULL

    filter <- apply(tab1,2,function(x) sum(x>0)/length(x)>input$min.pc.feature)
    tab1 <- tab1[,filter]
    tab1 <- tab1[,order(apply(tab1,2,function(x) sum(x>0)/length(x)),decreasing=T)]
    
    m1 <- match(id1,row.names(tab1))
    
    tab2 <- as.data.frame(data.chip[data.chip$ATAC.id%in%id2,])
    row.names(tab2) <- tab2$ATAC.id
    tab2$ATAC.id <- NULL

    tab2 <- tab2[,apply(tab2,2,function(x) sum(x>0)/length(x)>input$min.pc.feature)]
    tab2 <- tab2[,order(apply(tab2,2,function(x) sum(x>0)/length(x)))]
    m2 <- match(id2,row.names(tab2))

    tab.interact <- cbind(tab1[m1,],tab2[m2,])
    ids <-  paste(row.names(tab1)[m1],row.names(tab2)[m2],sep=":")
    tab.interact <- tab.interact[!duplicated(ids),]
    row.names(tab.interact) <- ids[!duplicated(ids)]
    annots <- data.frame(id1=row.names(tab1)[m1][!duplicated(ids)],id2=row.names(tab2)[m2][!duplicated(ids)],row.names=row.names(tab.interact))
    cut <- c(frag1=ncol(tab1),frag2=ncol(tab2))
    colnames(tab.interact) <-  unlist(sapply(strsplit(colnames(tab.interact),"[.][.]"),"[[",1))
    clust <-  pheatmap(t(tab.interact),
                       clustering_distance_cols = "binary",
                       clustering_method = "ward.D",silent=T)
    nb.clust <- input$nb.clust.heatmap.interact
    cluster=factor(cutree(clust$tree_col, k = nb.clust))
    annots$heatmap <- cluster[match(row.names(annots),names(cluster))]
    annots <- annots[match(names(cluster),row.names(annots)),]
    #annots$position <- 1:nrow(annots)
    
    list(tab.interact=tab.interact,annots=annots,cut=cut)
  })
  
  output$interacting.heatmap <- renderPlot({
    homer.tmp <- prepare.shaps()$homer.tmp
    homer.tmp$db.clust <- as.vector(homer.tmp$db.clust)
    homer.tmp$db.clust[!is.na(homer.tmp$cl1.split)] <- homer.tmp$cl1.split[!is.na(homer.tmp$cl1.split)] 
    data <- interacting.heatmap()
    annots <- data$annots
    data.interaction <- interacting.tSNE()
    annots$cl1 <- data.interaction$cluster1[match(annots$id1,data.interaction$ATAC.id1)]
    annots$cl2 <- data.interaction$cluster2[match(annots$id2,data.interaction$ATAC.id2)]
    annots$Feature <- data.interaction$class[match(annots$id1,data.interaction$ATAC.id1)]
    annots$frag1 <- peaks.ATAC.frag$frag.type[match(annots$id1,peaks.ATAC.frag$id)]
    annots$frag2 <- peaks.ATAC.frag$frag.type[match(annots$id2,peaks.ATAC.frag$id)]
    annots$status.1 <- peaks.ATAC.frag$status[match(annots$id1,peaks.ATAC.frag$id)]
    annots$status.2 <- peaks.ATAC.frag$status[match(annots$id2,peaks.ATAC.frag$id)]
    annots$genom.1 <- peaks.ATAC.frag$genomic[match(annots$id1,peaks.ATAC.frag$id)]
    annots$genom.2 <- peaks.ATAC.frag$genomic[match(annots$id2,peaks.ATAC.frag$id)]
    
    color.clust <- c(rep(c(brewer.pal(8, "Set1"),brewer.pal(8, "Set2"),brewer.pal(8, "Accent"),brewer.pal(8, "Spectral"),brewer.pal(8, "Pastel1"),brewer.pal(10, "BrBG")),3))
    names(color.clust) <- c(unique(homer.tmp$db.clust),"not.in.tSNE")
    color.clust <- color.clust[!is.na(names(color.clust))]
    color.clust <- color.clust[names(color.clust)%in%c(annots$cl1,annots$cl2)]
    pheatmap(t(data$tab.interact),show_rownames = T,
             clustering_method = "ward.D",
             annotation_col = annots[,rev(c("heatmap","Feature","cl1","cl2","frag1","frag2","genom.1","genom.2","status.1","status.2"))],
             color=colorRampPalette(colours()[c(235,34,35,24)])(100),
             border_color = NA,clustering_distance_cols = "binary",
             show_colnames=F,cluster_rows = F,gaps_row = data$cut[1],
             annotation_colors = list(cl1=color.clust[names(color.clust)%in%annots$cl1],cl2=color.clust[names(color.clust)%in%annots$cl2],genom.1=c(intergenic="grey",tss=colours()[614],prom=colours()[76],gene.body=colours()[613]),genom.2=c(intergenic="grey",tss=colours()[614],prom=colours()[76],gene.body=colours()[613]),frag1=c(o="black",b="white"),frag2=c(o="black",b="white"),Feature=c(gain_up="red",gain_down=colours()[53],loss_up=colours()[96],loss_down=colours()[124]),status.1=c(closing=colours()[131],opening="red",stable="grey"),status.2=c(closing=colours()[131],opening="red",stable="grey")))
  })
  
  output$select.on.heatmap.interaction.cluster <- renderUI({
    if(input$select.on.heatmap.interaction){
      data <- interacting.heatmap()
      annots <- data$annots
      checkboxGroupInput("heatmap.interaction.cluster","",choices = c("all",unique(annots$heatmap)),selected = "all",inline = T)
    }
  })
  
  select.on.frag <- reactive({
    data.interaction <- interacting.tSNE()
    #ids <- paste0(data.interaction$cluster1,"_",data.interaction$ATAC.id1,"_",data.interaction$fragment.id1,":",data.interaction$cluster2,"_",data.interaction$ATAC.id2,"_",data.interaction$fragment.id2)
    data.interaction$ATAC.id2[is.na(data.interaction$ATAC.id2)] <- "no.ATAC"
    id1 <- data.interaction$ATAC.id1
    id2 <- data.interaction$ATAC.id2
    ids <- paste(id1,id2,sep=":")
    if(input$select.on.heatmap.interaction){
      data <- interacting.heatmap()
      annots <- data$annots
      tab.interact <- data$tab.interact
      cut.cl <- data$cut
      if (input$heatmap.interaction.cluster[1]!="all"){
       annots <- annots[annots$heatmap%in%input$heatmap.interaction.cluster,]
       tab.interact <- tab.interact[row.names(tab.interact)%in%row.names(annots),]
      }
      frag1 <- tab.interact[,1:cut.cl[1]]
      frag1 <- colnames(frag1)[apply(frag1,2,function(x) sum(x>0)/length(x)>input$min.pc.feature)]
      frag2 <- tab.interact[,(cut.cl[1]+1):ncol(tab.interact)]
      frag2 <- c(colnames(frag2)[apply(frag2,2,function(x) sum(x>0)/length(x)>input$min.pc.feature)])
    }
    else {
      homer.num$Interaction_Number <- NULL
      homer.num <- merge(external.chip.ATAC,homer.num,by="ATAC.id")
      colnames(homer.num) <- unlist(sapply(strsplit(colnames(homer.num),"[.][.]"),"[[",1))
      
      frag1 <- homer.num[homer.num$ATAC.id%in%id1,]
      frag1 <- colnames(frag1)[apply(frag1,2,function(x) sum(x>0)/length(x)>input$min.pc.feature)][-1]
      frag2 <- homer.num[homer.num$ATAC.id%in%id2,]
      frag2 <- c(colnames(frag2)[apply(frag2,2,function(x) sum(x>0)/length(x)>input$min.pc.feature)][-1],"no.ATAC")
    }
  list(frag1=frag1,frag2=frag2)

  })
  
  output$select.on.frag1 <- renderUI({
    features <- select.on.frag()$frag1
    checkboxGroupInput("frag1.selected","",choices = c("all",features),selected = "all")
  })
  
  output$select.on.frag2 <- renderUI({
    features <- select.on.frag()$frag2
    checkboxGroupInput("frag2.selected","",choices = c("all",features),selected = "all")
  })
  
  tab.to.display <- reactive({
    heatmap.data <- interacting.heatmap()$annots
    data.interaction <- interacting.tSNE()
    data.interaction$ATAC.id2[is.na(data.interaction$ATAC.id2)] <- "no.ATAC"
    id1 <- data.interaction$ATAC.id1
    id2 <- data.interaction$ATAC.id2
    ids <- paste(id1,id2,sep=":")
    
    frag1.f <- select.on.frag()$frag1
    frag2.f <- select.on.frag()$frag2
    if (input$frag1.selected[1]!="all"){
      frag1.f <- input$frag1.selected
    }
    if (input$frag2.selected[1]!="all"){
      frag2.f <- input$frag2.selected
    }
    homer.num$Interaction_Number <- NULL
    homer.num <- merge(external.chip.ATAC,homer.num,by="ATAC.id")
    colnames(homer.num) <- unlist(sapply(strsplit(colnames(homer.num),"[.][.]"),"[[",1))
    homer.num <- as.data.frame(homer.num)
    frag1 <- homer.num[match(id1,homer.num$ATAC.id),]
    frag2 <- homer.num[match(id2,homer.num$ATAC.id),]
    frag2$ATAC.id[is.na(frag2$ATAC.id)] <- "no.ATAC"
    if (length(frag1.f)==1){
      keep.frag1 <- frag1[,frag1.f]>0
    }
    else{
      keep.frag1 <- apply(frag1[,frag1.f],1,function(x) sum(x>0))
      keep.frag1 <- ifelse(keep.frag1==length(frag1.f),TRUE,ifelse(keep.frag1==0,FALSE,ifelse(input$OR.AND.frag1=="OR",TRUE,FALSE)))
    }
    
    if ("no.ATAC"%in%frag2.f){
      frag2.f <- frag2.f[frag2.f!="no.ATAC"]
      frag2[is.na(frag2)] <- 1
    }
    else{
      frag2[is.na(frag2)] <- 0
    }
    if (length(frag2.f)==1){
      if (frag2.f=="no.ATAC"){
        keep.frag2 <- keep.frag2$ATAC.id=="no.ATAC"
      }
      else{
        keep.frag2 <- frag2[,frag2.f]>0
      }
    }
    else{
      keep.frag2 <- apply(frag2[,frag2.f],1,function(x) sum(x>0))
      keep.frag2 <- ifelse(keep.frag2==length(frag2.f),TRUE,ifelse(keep.frag2==0,FALSE,ifelse(input$OR.AND.frag2=="OR",TRUE,FALSE)))
    }
    data.interaction$interact.heatmap <- heatmap.data$heatmap[match(ids,row.names(heatmap.data))]
    
    data.interaction <- data.interaction[keep.frag1 & keep.frag2]
    if(input$select.on.heatmap.interaction){
      if (input$heatmap.interaction.cluster[1]!="all"){
        data.interaction <- data.interaction[interact.heatmap%in%input$heatmap.interaction.cluster]
      }
    }
    data.interaction
    
  })
  
  output$selected.interacting.fragments <- renderDataTable({
    tab.to.display()
  })
  
  interacting.tab <- eventReactive(input$go.save.interacting.tab,{
    withProgress(message = 'Dowloading files ...', {
      tab <- tab.to.display()
      write.table(tab,file.path(data.dir,"tables",paste0(input$prefix.save.interacting.tab,".tab")),sep="\t",row.names=F,quote=F)
    })
    text <- paste("File",paste0(input$prefix.save.interacting.tab,".tab"),"has been written in the table folder.")
  })
  
  output$download.interacting.done <- renderText({
    interacting.tab()
  })
  
  tab.expression.interacting <- reactive({
    genes <- interacting.tSNE()
    genes <- unique(c(genes$fragment.id1,genes$fragment.id2))
    genes.all <- genes[!grepl("^chr|^nb.",genes)]
    genes <- tab.to.display()
    genes <- unique(c(genes$fragment.id1,genes$fragment.id2))
    genes.selected <- genes[!grepl("^chr|^nb.",genes)]
    RNAseq <- RNAseq[ID%in%genes.all]
    
    time.course <- time.course[,list(ID=SYMBOL,T0.T20,T20.T1,T1.T2 ,T2.T4,T4.T24)]
    expr <- merge(RNAseq[,list(ID,selected=ifelse(ID%in%genes.selected,T,F),RNAseq=status)],time.course,by="ID",all.x=T)
    expr
  })
  
  output$tab.expression.interacting <- renderDataTable({
    tab.expression.interacting()
  })
  
  interacting.tab.expr <- eventReactive(input$go.save.interacting.expr.tab,{
    withProgress(message = 'Dowloading files ...', {
      tab <- tab.expression.interacting()
      write.table(tab,file.path(data.dir,"tables",paste0(input$prefix.save.interactingexpr.tab,".tab")),sep="\t",row.names=F,quote=F)
    })
    text <- paste("File",paste0(input$prefix.save.interactingexpr.tab,".tab"),"has been written in the table folder.")
  })
  output$download.interacting.expr.done <- renderText({
    interacting.tab.expr()
  })
  
})

