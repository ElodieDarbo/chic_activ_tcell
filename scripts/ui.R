library(shiny)
library(shinythemes)
library(shinycssloaders)
library(plotly)
library(DT)
library(shinybusy)

row <- function(...) {
  tags$div(class="row", ...)
}

col <- function(width, ...) {
  tags$div(class=paste0("span", width), ...)
}


shinyUI(
  fluidPage(navbarPage(theme = shinytheme("journal"),
                       add_busy_spinner(spin = "fading-circle"),
                       #titlePanel("Cell cycle entry during T-cell activation: a genomic perspective"),
                       tabPanel("Accueil",
                                
                                h1("Cell cycle entry during T-cell activation: a genomic perspective"),
                                verbatimTextOutput("acceuilText"),
                                h3("General"),
                                p("This application aims at help in exploring gene expression, chromatine accessibility and structure during the activation of T-cells."),
                                p("A spinning weel appears in the top right corner when the application is busy computing, loading, displaying data."),
                                p("Sometimes before some actions are made, errors are displayed in bold red because of missing data or something similar. Do not worry, it is normal and once the action is done everything goes to normal."),
                                p("I have made a LOT of tests but I may have missed things, in this case, an error is written in bold red. Usually, the app does not crash and an other action can be done. If nothing work anymore ... just write to me describing exactly what you did and a screeshot showing the error."),
                                h3("Genomic visualisation"),
                                p("In this tab you will be able to look at the interaction data per locus, bait symbol or identifier. 
                                The key word is 'patience', indeed many data are loaded which makes more or less heavy graphs.
                                The panels containing the interactions are interactive but not the ones with ChIP-seq data."),
                                h3("t-SNE clustering on SHAPs"),
                                #p("There are two tabs in this section. 't-SNE clustering on SHAPs' aims at exploring t-SNE results made from machine learning analysis of TFBS 
                                #mapping on ATAC peaks associated with differential interactions which involve differentially expressed genes.
                                #'ATAC + CHiC status vs logFC RNAseq' is informative and is inspired from Burren et al. Figure 3c."),
                                p("This section aims at exploring t-SNE results made from machine learning analysis of TFBS 
                                mapping on ATAC peaks associated with differential interactions which involve differentially expressed genes."),
                                #h4("t-SNE clustering on SHAPs"),
                                #"This section is a  bit complex at first seen.",
                                h4("1. The first thing to do is to choose a feature (presence/absence of a given TFBS or TF, ATAC peak annotation)."),
                                p("By default 'db.clust' is selected and represents the t-SNE clusters.
                                A check box allows to display the cluster ID, it needs to be checked before clicking on 'Show features'. It can be selected at any time but a re-click is needed to take it into account."),
                                p("At clicking on 'Show features', the t-SNE clustering with coloured peaks is diplayed. It has to be clicked each time you change of feature."),
                                h4("2. You are interested in a particular set of peaks:"),
                                h5("Select on the graph by click and drag"),
                                p("- This operation is very sensitive, any click on the graph is taken into account."),
                                p("- CAUTION ! If you select a new feature to display, be sure no area is selected in the plot. 
                                If one is selected, just click once in the plot outside of the area."),
                                p("- Bellow the graph, the number of selected peaks is displayed"), 
                                p("Be aware that all the peaks in the area selected are taken into account. 
                                Even the ones not allocated to any cluster and labelled '0'"),
                                h5("Focus on peaks in the area"),
                                p("Above the graph two check boxes are available: 'Focus on cluster' and 'Focus on feature'"),
                                p("- By checking the box 'Focus on cluster' you can select fragments from the most represented cluster in selected area."),
                                p("- By checking the box 'Focus on feature' you can select fragments according to the chosen feature in selected area.
                                In fact, clicking this box will display available value(s) for the displayed feature. 
                                If numeric, a slider appears (often 0 to 1), if categorical, all possible values are displayed and more than one can be checked."),
                                p("- Both Boxes can be checked and the peaks are automatically selected, you can see the number of peaks changing bellow the graph."),
                                p("- You can select the whole graph, focusing or not on features, both are very interesting."),
                                h4("3. Analyse selected peaks."),
                                p("Once you are happy with you selection, you can use the different tabs."),
                                p("Note that the last tab 'Heatmap TFBS/chip frequency clusters' is independent of the selected peaks and shows the complete set."),
                                h5("TF(BS)s frequencies"),
                                p("In this section appear 2 graphs. Left: the frequency of peaks with given TF(BS)s in a decreasing order. Right: heatmap showing the number of occurence of TF(BS)s in selected peaks with some annotations."),
                                p("- You have the possiblity to choose the TF(BS)s according to their frequency in the set of peaks using the slider 'Select x% of peaks'. This selection is applied to all subsections. By default, it is set at 20%."),
                                p("- You can split the peaks according to the clustering showed in the heatmap. It will label the peaks in tables of the other subsections."),
                                h5("Summary table for peaks"),
                                p("Details on selected ATAC peaks are displayed through a dynamic table where elements can be searched and columns can be ordered."),
                                p("You can download this table giving a name to the file in the text box and clicking 'Dowload table', an extention '.tab' is added to the file name by default."),
                                h5("Summary plots for peaks"),
                                p("Several plots show the distribution of peak annotations. If clusters are asked from the heatmap in 'TF(BS)s frequencies' tab, the plots are splitted accordingly."),
                                h5("Interacting fragments"),
                                p("- Details on selected ATAC peaks and their interacting fragments (with or whitout ATAC peaks) are displayed through a dynamic table where 
                                  elements can be searched and columns can be ordered. The heatmap.info column refers to the heatmap clusters from 'TF(BS)s frequencies'. Cluster1 and 2 columns refer to t-SNE clustering if applicable."),
                                p("- Distribution of ATAC peaks in interacting fragments in t-SNE clusters."),
                                p("- Heatmap showing occurences of TF(BS)s in interacting ATAC peaks. The rows of the heatmap are separated, on the top the TF(BS)s from the selected peaks (labelled in the annotations with suffix '1', e.g. cl1, genom.1 etc.), bottom part concerns the interacting ATAC peaks (labelled in the annotations with suffix '2', e.g. cl2, genom.2 etc.). Only TF(BS)s occurring in at least x% of the peaks ('TF(BS)s frequencies')."),
                                p("As in the 'TF(BS)s frequencies' subsection, you can choose to split the interactions in groups according to the clustering showed in the heatmap."),
                                p("- You can select particular interactions according to the TF(BS)s present in the ATAC peaks of each fragment ('no.ATAC' allows to keep interactions with fragments without ATAC peaks)."),
                                p("Checking the box 'Select on heatmap' will select ATAC peaks and their interaction from the heatmap, meaning interaction with fragment without ATAC peaks are discarded."),
                                p("By default all the interactions from the heatmap are selected. You can select one or more particular clusters from the heatmap by checking the box(es) that appeared in the right. You may deselect the box 'all' if choosing cluster(s)."),
                                p("- TF(BS)s present in at least x% of the ATAC peaks are displayed for both fragments. The OR (default) means that the ATAC peaks carry at least one of the selected TF(BS)s, AND means ATAC peaks have to carry all selected ones."),
                                p("By default, all TF(BS)s are selected (checkbox 'all' that has to be deselected if choosing to select)."),
                                p("- The result of the selection(s) is displayed with a table (similar to the one at the top of the section) with an added column 'interact.heatmap' corresponding to the clustering from the heatmap."),
                                p("You can download this table giving a name to the file in the text box and clicking 'Dowload table', an extention '.tab' is added to the file name by default."),
                                p("Some details about genes involved in the interactions are displayed. The column selected indicates if the gene is part of the final selection."),
                                p("You can download this table giving a name to the file in the text box and clicking 'Dowload table', an extention '.tab' is added to the file name by default."),
                                h3("Misc"),
                                p("In this section you have access to the data (interactions, fragments, ATAC peaks, gene expression) through interactive tables in which you can easily search for genes, fragments etc."),
                                p("You can dowload this tables, their names are predefined."),
                                h3("Enjoy !")
                                
                                
                       ),
                       tabPanel("Genomic visualisation",
                                h5("Computing and displaying the data can be a bit slow, the larger the area, the slower. Over 1 Mb, it is quite slow."),
                                textInput("coords.TADs","Enter genomic coordinates (chr:start-end), a gene symbol (capital letter) or a fragment identifier",value = "MIR155HG",width = 500),
                                h5("Please check your query before launching."),
                                actionButton("ckeck.exists","Check"),
                                textOutput("region.exists"),
                                selectInput("arc.plot","Select interaction type to plot",choices=c("differential","all")),
                                h5("NB1: if choosing differential interaction but no exists, stable ones are displayed."),
                                h5("NB2: if querying a gene or a fragment ID only interactions for this fragment are shown. If you want to display all interactions in the region, copy/paste the coordinates."),
                                actionButton("go.CHiC.plot.TAD","Apply"),
                                verbatimTextOutput("CHiC.TAD.plot.coords"),
                                plotOutput("chipSeq",width=1000,height = 500),
                                plotOutput("CHiC.TAD.plot",width=1000,
                                           brush = brushOpts(id = "plot2_zoom",resetOnNew = TRUE)),
                                h5("The plot showing interaction is zoomable. You can click and drag around the area you are interested in, it will zoom in the horizontal scale. If the displayed area is large, then the operation can take a few seconds. The larger, the slower. Caution! It is very sensitive so just a click has an impact."),
                                verbatimTextOutput(("CHiC.TAD.plot.coords.zoom")),
                                plotOutput("zoomin.chipSeq",width=1000,height = 500),
                                plotOutput("zoomin.CHiC",width=1000,click ="plot_hover"),
                                splitLayout(h5("Zoom in or out"),h5("Go left or right"),cellWidths = 300),
                                splitLayout(actionButton("zoomout","out",width = 50),actionButton("zoomin","in",width = 50),actionButton("Go.left","<<",width = 50),actionButton("Go.right",">>",width = 50),cellWidths = 145),
                                splitLayout(radioButtons("kb.zoom","Choose size of the zoom (in Kb)",choices = c(1,5,10,50,100,200),selected = 1,inline = T),radioButtons("kb.move","Choose size of the move (in Kb)",choices = c(1,5,10,50,100,200),selected = 100,inline = T),cellWidths = 300),
                                h5("You can click on the features (ATAC peaks, genes and fragments). It will show some more details in the table bellow."),
                                dataTableOutput("hover.info")
                       ),
                       tabPanel("t-SNE clustering on SHAPs",
                                #tabsetPanel(
                                  
                                  #tabPanel("t-SNE clustering on SHAPs",
                                           
                                             #tabPanel("Coloring graph", 
                                                      uiOutput("TFs.tsne.list"),
                                                      checkboxInput("show.db.clust","Show cluster number"),
                                                      actionButton("go.tSNE.shaps","Show features"),
                                                      h5("To analyse a particular subset of peaks, select an area on the plot and details on selected peaks will appear bellow. By default all peaks in the area are selected."),
                                                      h5("Showing the cluster number is not automatic, you need to check the box and re-click on 'Show features'."),
                                                      h5("CAUTION ! If you select a new feature to display, be sure no area is selected in the plot. If one is selected, just click once in the plot outside of the area."),
                                                      
                                                      h5("---"),
                                                      h5("By checking the box 'Focus on cluster' you can select fragments from the most represented cluster in selected area."),
                                                      h5("By checking the box 'Focus on feature' you can select fragments according to the chosen feature in selected area."),
                                                      splitLayout(checkboxInput("show.all.points.in.area","Focus on cluster"),
                                                       checkboxInput("show.feature.in.area","Focus on feature")),
                                                      uiOutput("select.features"),
                                                      h5("---"),
                                                      h5("NB: You can select the whole picture."),
                                                      plotOutput("tSNE.shaps.icis",width=800,height=750,
                                                                 brush = brushOpts(id = "plot_brush")),
                                                      
                                                      textOutput("nb.frag.selected"),
                                                      h5(""),

                                                      
                                             #),
                                           tabsetPanel(
                                             tabPanel("TF(BS)s frequencies",
                                                      h5("TF(BS)s present in more than x% of the fragments"),
                                                      sliderInput("min.pc.feature",label = "Select x% of peaks",min = 0,max = 1,value = 0.2,step = 0.1),
                                                      h5("If you want to split the peaks according to the heatmap, choose the number of clusters. This information will appear in the summary table."),
                                                      sliderInput("nb.clust.heatmap",label = "Number of clusters",min = 1,max = 10,value = 1,step = 1),
                                                      splitLayout(plotOutput("feature.shaps.clust",width=500,height=1000),
                                                      
                                                      plotOutput("heatmap.shaps.clust.plot",height=1000,width=500))
                                                      
                                             ),
                                             tabPanel("Summary table for peaks",
                                                      dataTableOutput("brush_info"),
                                                      h5("Prefix for the name of the table to dowload. A column for each TF(BS) present in the heatmap is added."),
                                                      textInput("prefix.save.shaps.tsne", "", value = ""),
                                                      actionButton("go.save.shaps.tsne","Dowload table"),
                                                      textOutput("download.done")
                                        
                                             ),
                                             tabPanel("Summary plots for peaks",
                                                      plotOutput("tab.shaps.clust.plot",height=800),
                                                      plotOutput("tab.shaps.clust.plot2",height=400),
                                                      plotOutput("tab.shaps.clust.plot3",height=800)
                                                      
                                                      
                                             ),
                                             tabPanel("Interacting fragments",
                                                      dataTableOutput("tab.interacting.tSNE"),
                                                      plotlyOutput("all.shaps.clust"),
                                                      sliderInput("nb.clust.heatmap.interact",label = "Number of clusters",min = 1,max = 10,value = 1,step = 1),
                                                      plotOutput("interacting.heatmap",height=1000),
                                                      h5("Select interactions according to feature on fragments"),
                                                      splitLayout(checkboxInput("select.on.heatmap.interaction","Select on heatmap"),uiOutput("select.on.heatmap.interaction.cluster")),
                                                      splitLayout(radioButtons("OR.AND.frag1","Fragment with selected ATAC peaks",choices = c("OR","AND"),selected = "OR"),radioButtons("OR.AND.frag2","Interacting fragments",choices = c("OR","AND"),selected = "OR")),
                                                      splitLayout(uiOutput("select.on.frag1"),uiOutput("select.on.frag2")),
                                                      dataTableOutput("selected.interacting.fragments"),
                                                      h5("Prefix for the name of the table to dowload. "),
                                                      textInput("prefix.save.interacting.tab", "", value = ""),
                                                      actionButton("go.save.interacting.tab","Dowload table"),
                                                      textOutput("download.interacting.done"),
                                                      h5("Dowdload the expression table for further functional analyses (e.g. https://maayanlab.cloud/Enrichr/ , http://biit.cs.ut.ee/gprofiler/gost)"),
                                                      dataTableOutput("tab.expression.interacting"),
                                                      textInput("prefix.save.interactingexpr.tab", "", value = ""),
                                                      actionButton("go.save.interacting.expr.tab","Dowload table"),
                                                      textOutput("download.interacting.expr.done")
                                                      

                                           ),
                                           tabPanel("Heatmap TFBS/chip frequency clusters",
                                                    plotOutput("heatmap.freq.TFBS",width=1200,height=600),
                                                    plotOutput("heatmap.freq.chip",height=600),
                                                      )
                                           #tabPanel("Cluster network",
                                                    
                                          # ),
                                           #tabPanel("Cluster mapping other datasets",
                                                    
                                           #)
                                           )
                                           
                                   #) ,
                                   
                                  # tabPanel("ATAC + CHiC status vs logFC RNAseq",
                                  #          plotOutput("ATAC.bait.RNA",width=1000,heigh=800)
                                  # )
                                #)
                       ),
                       tabPanel("Misc",
                                tabsetPanel(
                                  tabPanel("Interactions",
                                           checkboxInput("only.diff","Select differential"),
                                           dataTableOutput("tab.interaction"),
                                           actionButton("save.tab.int","Download"),
                                           textOutput("write.tab.interaction")
                                  ),
                                  tabPanel("Fragments",
                                           checkboxInput("only.diff.frag","Select fragments with differential interaction(s)"),
                                           dataTableOutput("tab.fragments"),
                                           actionButton("save.tab.frag","Download"),
                                           textOutput("write.tab.fragments")
                                  ),
                                  tabPanel("ATAC peaks",
                                           dataTableOutput("tab.ATAC"),
                                           actionButton("save.tab.ATAC","Download"),
                                           textOutput("write.tab.ATAC")
                                  ),
                                  tabPanel("Expression",
                                           dataTableOutput("tab.expression"),
                                           actionButton("save.tab.expr","Download"),
                                           textOutput("write.tab.expression")
                                  )
                                )
                                
                       )
          )      
  )
  
)