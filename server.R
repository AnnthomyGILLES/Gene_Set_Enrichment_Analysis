#
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

# Importer library
library(shiny)
# library(shinydashboard)
library(clusterProfiler)
library(AnnotationDbi)
library(GO.db)
library(topGO)
library(clusterProfiler)
library(ReactomePA)
# library(biomaRt) 
library(DT)
library(shinythemes)
# library(shinyjs)
library(shinyBS)
library(shinycssloaders)
library(DOSE)
library(pathview)
library(png)

#Load Organism DB
library("org.Hs.eg.db")

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
  
  # Set constents ==============================================================================
  keggDirectory="./KeggPathways"
  keggResults="./KeggResults"
  
  # ############################################################################################
  #                                       General Functions
  # ############################################################################################
  
  # Data Selection =============================================================================
  
  # Gene Id column selector--------------------------------------------------------------------- 
  IdSelector<-function(data, inputId) {
    if(is.null(data) ) {
      return("No file uploaded.")
    } else {
      # Id Selection
      return(selectInput(inputId = inputId, label = "Select feature name :",
                  choices = colnames(data)))
    }
  }
  
  # Ordering/selecting criterion selector -------------------------------------------------------
  CritSelector<-function(data, inputId) {
    if(is.null(data) ) {
      return("No file uploaded.")
    } else {
      # Id Selection
      return(selectInput(inputId = inputId, label = "Select criterion :",
                         choices = colnames(data)))
    }
  }
  
  # Select direction --------------------------------------------------------------------------------------
  selectDir<-function(test, data, inputId){
    if(test=="SEA"){
      if(! is.null(data)) {
        return(checkboxInput(inputId= inputId, label = "Filter gene by upper criterion :", value = FALSE))
      } else { return(NULL)}
    }
  }
  
  # CutOff Selection --------------------------------------------------------------------------------------
  cutOffSelector<-function(test, data, colName, inputId ) {
    if(test=="SEA"){
      if(! is.null(data)) {
        # Exctract criterion column
        vec<-as.numeric(data[,grep( colName, colnames(data))])
        # Clean vector
        if((sum(is.na(vec))>0)) vec<-vec[-which(is.na(vec))]
        if(sum(abs(vec)==Inf)>0) vec<-vec[-which(abs(vec)==Inf)]
        # Construct selector
        return(sliderInput(inputId = inputId, label = "Select cutoff for gene set selection :",
                    min = min(vec), max = max(vec), value = median(vec), 
                    step = abs(max(vec)-min(vec))/1000, dragRange =FALSE))
        
      } else { return(NULL)}
    }
  }
  
  # GeneSet selection ---------------------------------------------------------------------------------------
  geneSubset<-function(test, data, id, criterion, value, upper=TRUE) {
    
    #get indices
    critCol<-grep(criterion, colnames(data))
    idCol<-grep(id, colnames(data))
   
    if(test == "SEA") {
      if(upper) set<-data[which(data[, critCol]>value), idCol]
      else set<-data[which(data[, critCol]<=value), idCol]
      return(as.character(set))
    } else {
      set<-data()[,critCol]
      names(set)<-as.character(data[,idCol])
      set<-set[order(set, decreasing=TRUE)]
      return(set)
    }
  }
  
  # Tests Settings ============================================================================
  
  # Set Max P value ---------------------------------------------------------------------------
  pValCutoff<-function(inputId, min=0.001, max=0.10, value=0.05, step=0.001){
    return(sliderInput(inputId = inputId , label = "Choose P-value cutoff:",
                min = min, max = max, value = value, step = step,
                dragRange =FALSE))
  }
  
  # Set Max Q value/Number of permutations ----------------------------------------------------
  testOption2<-function(test, inputIdQCutoff, inputIdNPerm, 
                        minQVal=0.05, maxQVal = 0.30, defaultQVal=0.10, QValStep = 0.05,
                        minNPerm= 10, maxNPerm=1e6, defaultNPerm=500, NPermStep=10) {
    if(test =="SEA") {
      # Q-value cutoff
      return(sliderInput(inputId = inputIdQCutoff, label = "Choose Q-value cutoff:",
                        min = minQVal, max = maxQVal, value =defaultQVal , step = QValStep,
                        dragRange =FALSE))
    } else {
      # Number of permutations
      return(sliderInput(inputId = inputIdNPerm, label = "Choose number of permutations :",
                        min =minNPerm, max = maxNPerm, value = defaultNPerm, step = NPermStep,
                        dragRange =FALSE))
    }
  } 
    
  # Set gene set sizes ------------------------------------------------------------------------ 
  selectSize<-function(inputId, min=10, max=800, step=10, value=c(100, 500)) {
    # check value
    if (length(value) != 2) { warning("2 values must be set")}
    # min and max Gene Set size
    return(sliderInput(inputId = inputId, label = "Choose set min and max size :",
                min = min, max = max, step = step, value = value,
                dragRange =TRUE))
  }  
  
  # Set correction ----------------------------------------------------------------------------
  setCorrection<-function(inputId, choices=list( "bonferroni", "BH", "none"), selected="BH") {
    return(selectInput(inputId = inputId, label = "Select correction to run :",
                choices = choices, selected = selected))
  }
  
  
  # Plotting functions ========================================================================
  # Dot Plot for gse --------------------------------------------------------------------------
  dotplot.gse<-function(gseResult, title) {
    data<-slot(gseResult, "result")
    # order data
    data<-data[order(data$p.adjust, decreasing=FALSE),]
    # Set color ramp
    pal<-colorRampPalette(c("red", "blue"))
    col<-pal(length(unique(data$p.adjust)))[sapply(data$p.adjust, function(p){which(unique(data$p.adjust==p))})]
    # Select top results
    data<-data[1:min(10, nrow(data)),]
    # Set plotting parameters for legends
    par(mar=c(5,max(sapply(data$ID, nchar))/1.5, 4,6), xpd=TRUE)
    # Plot points
    plot(data$enrichmentScore, 1:nrow(data)+1, col=col, pch=19, 
         yaxt="n", ylab="", xlab="Enrichment Score", 
         cex=log10(data$setSize)+1, main= title)
    axis(2, at=1:nrow(data)+1, label=data$ID, las=2)
    # Add legends
    legend("topright", inset=c(-0.25, -0.04), 
           legend =round(unique(c(min(data$p.adjust), median(data$p.adjust), max(data$p.adjust))),3), 
           pch = 19, 
           col = pal(length(unique(data$p.adjust)))[sapply(unique(c(min(data$p.adjust),median(data$p.adjust), max(data$p.adjust))), function(level){
              result<-which(unique(data$p.adjust)==level)
              if(length(result)==0)return(0)
              else return(result[1])
           })],
           title="P value", bty="n")
    sizes<-unique(data$setSize)
    sizes<-sizes[order(sizes)]
    legend("topright", inset=c(-0.25, 0.4), legend=sizes, pch=19, 
           col="gray", pt.cex=log10(sizes)+1, bty="n", title="Set Size")
  }
  
  
  # Integrated dot plot funtion ---------------------------------------------------------------
  dotplotGO<-function(data, test, onto) {
    if(nrow(slot(data, "result"))==0) return(warning("No enriched term"))
    
    if(test== "SEA") dotplot(data)
    else dotplot.gse(data, paste0("Top 10 enriched", onto, "terms"))
  }
  
  # Integrated GO plot function ---------------------------------------------------------------
  GOplot<-function(envResult) {
    # extract  data
    data<-envResult
    # extract table
    table<-slot(data, "result")
    # Remove NAs
    i<-grep("NA", table$ID)
    if(length(i)>0) table<-table[-i,]
    # Check number of rows
    if(nrow(table)==0) return(warning("No enriched term"))
    # Replace table
    slot(data, "result")<-table
    # Plot
    plotGOgraph(data, firstSigNodes = min(10, nrow(table)))
  }
  
  # Treat input data ==========================================================================
  
  # Import Data--------------------------------------------------------------------------------
  data<-reactive({
    file<-input$inFile
    if(! is.null(file)) {
      read.table(file$datapath, header=input$header, sep=input$sep,quote = input$quote, stringsAsFactors=FALSE)
    }else {
      return(NULL)
    }
  })
  
  # Shiny Dashobard----------------------------------------------------------------------------
  output$eventTable <- renderDataTable({
    data()
    
  })
  
  
  # ############################################################################################
  #                                        GENE ONTOLOGY
  # ############################################################################################ 
  
  # Save test choice----------------------------------------------------------------------------
  test<-reactive(input$test)
  
  # Gene Id Selector--------------------------------------------------------------------------
  output$IdSelector<-renderUI({
    IdSelector(data(), "geneId")
  })
  
  
  # Criterion Selector--------------------------------------------------------------------------
  output$OrderSelector<-renderUI({
    CritSelector(data(), "criterion")
  })

  
  # Test Selector--------------------------------------------------------------------------------
  output$upper<-renderUI({
    selectDir(test(), data(), "upperSelect")
  })

  # Gene Set Cutoff Selector--------------------------------------------------------------------------------------
  output$setCutoff<-renderUI({
    cutOffSelector(test(), data(), input$criterion, "setSelect")
  })
  
  # Set up tests=================================================================================
  
  
  
  # Annotations options--------------------------------------------------------------------------
 
  # Choose to activate GO level filtering
  output$annOption2<-renderUI({
    # Activate GO level filtering
    checkboxInput(inputId= "filter", label = "Filter by GO level :", value = FALSE)
  })
  
  # Select GO level filtering option
  output$annOption3<-renderUI({
    if(! is.null(input$filter)) {
      if(input$filter) {
        # GO level selection
        numericInput(inputId = "level", label = "Select level to consider :",
                     min = 1, max = 12, step = 1, value = 3)
      }
    }
  })
  
  
  # Choose significativity threshold -----------------------------------------------------------
  output$testpCutoff<-renderUI({
    pValCutoff("pCutoff")
  })

  # Choose q value threshold -------------------------------------------------------------------
  output$testqCutoff<-renderUI({
    testOption2(test(), "qCutoff", "nPerm")
  })
  
  
  # Choose minimum and maximum size to consider a genes set -----------------------------------
  output$testSetSize<-renderUI({
    selectSize("setSize")
  })

  # Choose correction method for multiple testing ---------------------------------------------
  output$testCorrection<-renderUI({
    setCorrection("correction")
  })
  
  
  
  # Run test======================================================================================================================================================
  
  # Get geneSet --------------------------------------------------------------------------------------------------------------------------------------------------
  
  geneSet<-reactive({
    return(geneSubset(test(), data(), input$geneId, input$criterion, input$setSelect, input$upperSelect))
  })
  
  
  
  # Plot GO Set distributions=====================================================================================================================================

  # Calcualte GO distribution function ----------------------------------------------------------------------
  gogroup<-function(test, onto,level=2) {
    if(test =="SEA") {
      # Calculate whole Go Set distribution
      go <- groupGO(gene = geneSet(),
                    OrgDb = 'org.Hs.eg.db', #! create organism selector
                    ont= onto,
                    level=level,
                    keytype="ENSEMBL",
                    readable=TRUE)
      return(go)
    } else {
      genes<-c()
      genes<-names(geneSet())
      # else genes<-geneSet()
      go <- groupGO(gene = genes,
                    OrgDb = 'org.Hs.eg.db', #! create organism selector
                    ont= onto,
                    level=level,
                    keytype="ENSEMBL",
                    readable=TRUE)
      return(go)
    }
  }
  
  # Treat go distribution plotting ------------------------------------------------------------------
  distri <- eventReactive(input$plotWholeGO, {
      if(!input$filter){
        if (input$donum1) goBP<-gogroup(test(), "BP",2) else goBP<-NULL
        if (input$donum2) goCC<-gogroup(test(), "CC",2) else goCC<-NULL
        if (input$donum3) goMF<-gogroup(test(),"MF",2) else goMF<-NULL
      }else {
        if (input$donum1) goBP<-gogroup(test(), "BP",input$level) else goBP<-NULL
        if (input$donum2) goCC<-gogroup(test(), "CC",input$level) else goCC<-NULL
        if (input$donum3) goMF<-gogroup(test(), "MF",input$level) else goMF<-NULL
      }
      list(BP=goBP,CC=goCC,MF=goMF)
  })
 
  output$wholeGOSetBP = renderPlot({
    if(!is.null(distri()$BP) && input$donum1) {
        if(input$test=="SEA"){barplot(distri()$BP,  font.size = 6, drop=TRUE, showCategory=10, order=TRUE,title="GO BP Term Distribution in whole dataset")}
        else {barplot(distri()$BP,  font.size = 6, drop=TRUE, showCategory=10, order=TRUE,title="GO BP Term Distribution")}
    }
   })
  
  output$wholeGOSetCC = renderPlot({
    if(!is.null(distri()$CC) && input$donum2) {
      if(input$test=="SEA"){barplot(distri()$CC,  font.size = 6, drop=TRUE, showCategory=10, order=TRUE,title="GO CC Term Distribution in whole dataset")}
      else {barplot(distri()$CC,  font.size = 6, drop=TRUE, showCategory=10, order=TRUE,title="GO CC Term Distribution")}
    }
  })
  
  output$wholeGOSetMF = renderPlot({
    if(!is.null(distri()$MF) && input$donum3) {
      if(input$test=="SEA"){barplot(distri()$MF,  font.size = 6, drop=TRUE, showCategory=10, order=TRUE,title="GO MF Term Distribution in whole dataset")}
      else {barplot(distri()$MF,  font.size = 6, drop=TRUE, showCategory=10, order=TRUE,title="GO MF Term Distribution")}
    }
  })
  
  
  # output$downloadData <- downloadHandler({
  #   filename = "Shinyplot.png",
  #   content = function(filename) {
  #     png(filename,width = 1800, height = 900)
  #     plot(barplot(model()$BP,  font.size = 6, drop=TRUE, showCategory=10, order=TRUE,title="GO BP Term Distribution in whole dataset"))
  #     dev.off()
  #   }) 
  

  # Activate test ---------------------------------------------------------------------------------
  output$runTest<-renderUI({
    actionButton(inputId="runGO", label="Run Test :")
  })

  # Test function ---------------------------------------------------------------------------------
  goResult<-function(onto) {
    if(test() == "SEA") {
      data<-geneSet()
      return ( enrichGO(gene          = data,
                        #universe      = geneList,
                        keytype= "ENSEMBL",
                        OrgDb         = 'org.Hs.eg.db',
                        ont           = onto,
                        minGSSize    = input$setSize[1],
                        maxGSSize    = input$setSize[2],
                        pAdjustMethod = input$correction,
                        pvalueCutoff  = input$pCutoff,
                        qvalueCutoff  = input$qCutoff,
                        readable      = TRUE))
    } else {
      return(gseGO(geneList     = geneSet(),
                   keytype = "ENSEMBL",
                   OrgDb        = 'org.Hs.eg.db',
                   ont          = onto,
                   nPerm        = input$nPerm,
                   minGSSize    = input$setSize[1],
                   maxGSSize    = input$setSize[2],
                   pvalueCutoff = input$pCutoff,
                   pAdjustMethod = input$correction,
                   verbose      = FALSE)
             
      )
    }
  }

  # Gather Results --------------------------------------------------------------------------------
   enrichResult<-eventReactive(input$runGO, {
      if (input$donum1) goBP<-goResult("BP") else goBP<-NULL
      if (input$donum2) goCC<-goResult("CC") else goCC<-NULL
      if (input$donum3) goMF<-goResult("MF") else goMF<-NULL
      list(BP=goBP,CC=goCC,MF=goMF)
  })
  
   
  # Plot Dotplot  ---------------------------------------------------------------------------------
  output$dotplotBP<-renderPlot({
    if(!is.null(enrichResult()$BP) && input$donum1) dotplotGO(enrichResult()$BP, test(), "BP")
  })
  
  output$dotplotCC<-renderPlot({
    if(!is.null(enrichResult()$CC) && input$donum2) dotplotGO(enrichResult()$CC, test(), "CC")
  })
  
  output$dotplotMF<-renderPlot({
    if(!is.null(enrichResult()$BP) && input$donum3) dotplotGO(enrichResult()$MF, test(), "MF")
  })

  # Plot GOgraph -----------------------------------------------------------------------------------
  output$plotGOgraphBP<-renderPlot({
    if(!is.null(enrichResult()$BP) && input$donum1 && test()=="SEA") {
      GOplot(enrichResult()$BP)
    }
  })
  
  output$plotGOgraphCC<-renderPlot({
    if(!is.null(enrichResult()$CC) && input$donum2 && test()=="SEA") {
      GOplot(enrichResult()$CC)
    }
  })
  
  output$plotGOgraphMF<-renderPlot({
    if(!is.null(enrichResult()$MF) && input$donum3 && test()=="SEA") {
      GOplot(enrichResult()$MF)
    }
  })
  
  # Enrichmap
  # output$enrichMap<-renderPlot({
  #   if(! is.null(goResult())) {
  #     enrichMap(goResult(), n = 15)
  #   }
  # })
  


  
  # ############################################################################################
  #                                        PATHWAY
  # ############################################################################################ 
  
  # Import data ================================================================================
  
  # Save test option
  testPathway<-reactive(input$testPathway)
  
  # Select gene ID column ----------------------------------------------------------------------
  output$IdSelectorPathway<-renderUI({
    IdSelector(data(), "geneIdP")
  })
  
  # Select Ordering Criterion column ----------------------------------------------------------
  output$OrderSelectorPathway<-renderUI({
    CritSelector(data(), "criterionP")
  })
  
  # Set upper/lower selection -----------------------------------------------------------------
  output$upperPathway<-renderUI({
    selectDir(testPathway(), data(), "upperSelectP")
  })
  
 # Select cutoff --------------------------------------------------------------------------------
  output$setCutoffPathway<-renderUI({
    cutOffSelector(testPathway(), data(), input$criterionP, "setSelectP")
  })
  
  # Database selection --------------------------------------------------------------------------
  output$annOption1Pathway<-renderUI({
    # Database selection
    selectInput(inputId = "database", label = "Select pathway database :",
                choices = list( "KEGG", "DAVID", "REACTOME"))
  })
  
  # Tests Settings ===============================================================================
  
  # Choose significativity threshold -------------------------------------------------------------
  output$testpCutoffPathway<-renderUI({
    pValCutoff("pCutoffP")
  })
  
    # Choose q value threshold or nperm ----------------------------------------------------------
  output$testqCutoffPathway<-renderUI({
    testOption2(testPathway(), "qCutoffP", "nPermP")
  })
  
  # Choose minimum and maximum size to consider a genes set  -------------------------------------
  output$testSetSizePathway<-renderUI({
    selectSize("setSizeP", value=c(10,120))
  })
  
  # Select P-value correction --------------------------------------------------------------------
  output$testCorrectionPathway<-renderUI({
    setCorrection("correctionP")
  })
  
  # Run tests ====================================================================================
  # Run button -----------------------------------------------------------------------------------
  output$runTestPathway<-renderUI({
    actionButton(inputId="runPathway", label="Run Test :")
  })
  
  # Extract GeneSet ------------------------------------------------------------------------------
  geneSetPathway<-reactive({
    return(geneSubset(testPathway(), data(), input$geneIdP, input$criterionP, input$setSelectP, input$upperSelectP))
  })
  
 
  # Calculate enrichment Results ----------------------------------------------------------------
  PathwayResult<-eventReactive(input$runPathway,{
    if(testPathway() == "SEA") { 
      
      # Treat DAVID
      if(input$database=="DAVID") {
        if(length(geneSetPathway())>3000) set<-geneSetPathway()[1:3000]
        else set<-geneSetPathway()
        return(enrichDAVID(gene          = set,
                           #universe      = geneList,
                           idType= "ENTREZ_GENE_ID", # see for ID conversion
                           listTyp = "Gene",
                           annotation         = "KEGG_PATHWAY",
                           minGSSize    = input$setSizeP[1],
                           maxGSSize    = input$setSizeP[2],
                           pAdjustMethod = input$correctionP,
                           pvalueCutoff  = input$pCutoffP,
                           qvalueCutoff  = input$qCutoffP,
                           david.user = "annthomy.gilles@etu.univ-rouen.fr")
        )
      } else if(input$database=="KEGG") {
        kk <-enrichKEGG(gene          = geneSetPathway(),
                        #universe      = geneList,
                        keyType= "ncbi-geneid", # see for ID conversion
                        organism         = "hsa",
                        minGSSize    = input$setSizeP[1],
                        maxGSSize    = input$setSizeP[2],
                        pAdjustMethod = input$correctionP,
                        pvalueCutoff  = input$pCutoffP,
                        qvalueCutoff  = input$qCutoffP)
        return(kk)
      } else if(input$database=="REACTOME"){
        return ( enrichPathway(gene          = geneSetPathway(),
                               #universe      = geneList,
                               # keytype= "ENSEMBL",
                               organism         = "human",
                               minGSSize    = input$setSizeP[1],
                               maxGSSize    = input$setSizeP[2],
                               pAdjustMethod = input$correctionP,
                               pvalueCutoff  = input$pCutoffP,
                               qvalueCutoff  = input$qCutoffP,
                               readable      = TRUE)
        )
      } else { return(NULL)}
    }else {
      
      # Treat DAVID
      if(input$database=="DAVID") {
        warning("GSEA is not encoded")
        return(NULL)
      }else if(input$database =="REACTOME" ) {
        set<-geneSetPathway()
        names(set)<-as.character(names(set))
        return(gsePathway(geneList     = geneSetPathway(),
                          # keytype = "ENSEMBL",
                          organism        = "human",
                          nPerm        = input$nPermP,
                          minGSSize    = input$setSizeP[1],
                          maxGSSize    = input$setSizeP[2],
                          pvalueCutoff = input$pCutoffP,
                          pAdjustMethod = input$correctionP,
                          verbose      = FALSE))
      }else {
        return(gseKEGG(geneList     = geneSetPathway(),
                       keyType = "ncbi-geneid",
                       organism        = "hsa",
                       nPerm        = input$nPermP,
                       minGSSize    = input$setSizeP[1],
                       maxGSSize    = input$setSizeP[2],
                       pvalueCutoff = input$pCutoffP,
                       pAdjustMethod = input$correctionP,
                       verbose      = FALSE)
        )
      }
    }})
  
  # Graphical output =====================================================================
  # Dotplot ------------------------------------------------------------------------------
  output$dotplotPathway<-renderPlot({
    if(! is.null(PathwayResult())) {
      if(nrow(slot(PathwayResult(), "result"))==0) return(warning("No enriched pathway"))
      dotplot(PathwayResult())
    }
  })
  
  # Dotplot
  output$dotplotPathwayGSEA<-renderPlot({
    if(! is.null(PathwayResult())) {
      if(nrow(slot(PathwayResult(), "result"))==0) return(warning("No enriched pathway"))
      dotplot.gse(PathwayResult(), title = "Top 10 enriched pathways")
    }
  })
  
  # Barplot -------------------------------------------------------------------------------
  output$barplotPathway<-renderPlot({
    if(! is.null(PathwayResult())) {
      if(testPathway()== "SEA") {
        if(nrow(slot(PathwayResult(), "result"))==0) return(warning("No enriched pathway"))
        barplot(PathwayResult(), showCategory=10)
      }
    }
  })
  
  
  # Enrichmap ------------------------------------------------------------------------------
  output$enrichMapPathway<-renderPlot({
    if(! is.null(PathwayResult()) && testPathway()== "SEA") {
        if(nrow(slot(PathwayResult(), "result"))==0) return(warning("No enriched pathway"))
        enrichMap(PathwayResult(), n = 15)
      }
  })
  
  output$enrichMapPathwayGSEA<-renderPlot({
    if(! is.null(PathwayResult())&& testPathway()== "GSEA") {
        if(nrow(slot(PathwayResult(), "result"))==0) return(warning("No enriched pathway"))
        enrichMap(PathwayResult(), n = 15)
      }
  })
  
  
  # PathwayPlot --------------------------------------------------------------------------------
  # Get plot list
  pathwayList<-eventReactive(input$runPathway,
    # Check conditions
    if(! is.null(PathwayResult()) && input$database=="KEGG") {
      # Check results
      if(nrow(slot(PathwayResult(), "result"))==0) return(warning("No enriched pathway"))
      # Check ref dir
      if(! dir.exists (keggDirectory)) dir.create(keggDirectory)
      # Get enriched pathways
      idList<-slot(PathwayResult(), "result")$ID
      # Exctract gene info
      geneL<-geneSetPathway()
      if(testPathway()=="SEA") {
        geneL<-seq(1:length(geneSetPathway()))
        names(geneL)<-geneSetPathway()
      }
      print(geneL)
      print(idList)
      # Check result dir
      if(! dir.exists (keggResults)) dir.create(keggResults)
      # Plot
      pathview(gene.data=geneL, kegg.dir=keggDirectory, pathway.id=idList, species="hsa")
      # Get plot list
      plotList<-list.files(pattern="pathview")
      # Copy files
      file.copy(from=plotList, to=keggResults, overwrite=TRUE)
      # Let copy end
      Sys.sleep(0.2)
      # Remove files
      file.remove(plotList)
      # return pathways list
      return(idList)
    } 
  )
  
  # Give list to selector
  output$pathwayResultIds<-renderUI({
    if(! is.null(pathwayList()) && input$database=="KEGG" ) {
      return(selectInput(inputId = "pathviewID", label = "Select pathway to plot :",
                         choices = as.list(c(NULL, pathwayList()))))
      
    }
  })
  
  # Get reactive value
  pathviewID<-reactive(input$pathviewID)
  
  # Show plot
  output$pathwayPlot<-renderPlot({
    if(! is.null(pathviewID()) && input$database=="KEGG" && testPathway()=="SEA") {
      file<-paste0(keggResults, "/", pathviewID(), ".pathview.png" )
      if(file.exists(file)) {
        img<-readPNG(file)
        grid::grid.raster(img)
      }
    }
    
  })
  
  output$pathwayPlotGSEA<-renderPlot({
    if(! is.null(pathviewID()) && input$database=="KEGG" && testPathway()=="GSEA") {
      file<-paste0(keggResults, "/", pathviewID(), ".pathview.png" )
      if(file.exists(file)) {
        img<-readPNG(file)
        grid::grid.raster(img)
      }
    }
    
  })
  
  # ############################################################################################
  #                                        DOMAIN
  # ############################################################################################ 
  
  # Select data ================================================================================
  # Select gene ID column
  output$IdSelectorDomain<-renderUI({
    IdSelector(data(),"geneIdDomain")
  })
  
  # Select criterion column
  output$OrderSelectorDomain<-renderUI({
    CritSelector(data(), "criterionDomain")
  })
  
  # Subsetting direction
  output$upperDomain<-renderUI({
    selectDir("SEA", data(), "upperSelectDomain")
  })
  
  # 
  output$setCutoffDomain<-renderUI({
    cutOffSelector("SEA", data(), input$criterionDomain, "setSelectDomain")
  })
  
  # Set test options ===================================================================================
  # Set correction method
  output$testCorrectionDomain<-renderUI({
    setCorrection("correctionDomain", choices = list( "bonferroni", "BH"))
  })
  
  output$runTestDomain<-renderUI({
    actionButton(inputId="runDomain", label="Run Test :")
  })
  
  # Specific functions ==================================================================================
  mart<-useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  #mart<-NULL
  # Suppress emptyt answers -----------------------------------------------------------------------------
  clean.blankLines <- function(vector) #nettoyer les ids vides parfois renvoyés par BIOCONDUCTOR
  {
    return(vector[-which(vector=="")])
  }
  
  # Test function ---------------------------------------------------------------------------------------
  enrich.domains<-function(samples, domainIDs, universe, p.adj.method="BH") {
    
    #calculate p-values and nbhits
    results<-sapply(domainIDs, function(ID) {
      
      nbhits= sum(samples == ID)
      univhits = sum(universe == ID)
      univlength = length(universe)
      samplelength = length(samples)
      
      #enrichment test
      pval <-
        phyper(
          nbhits,#id nb hits in sample
          univhits, #id nb hits in universe
          univlength-univhits, #universe length - id hits in unbiverse
          samplelength, #sample length
          lower.tail = FALSE # test on rigth tail
        )
      
      return(c( nbhits=nbhits,
                expected.nbhits = univhits/univlength*samplelength, p.value=pval))
    })
    results<-data.frame(t(results))
    #adj p values
    results$p.adj <-p.adjust(results$p.value, method = p.adj.method, n = nrow(results)) 
    # order results
    results<-results[order(results$p.adj, decreasing=FALSE),]
    return(results)
  }
  
  
  # Doplot function --------------------------------------------------------------------------------------------------------
  dotplot.domains<-function(results, listGenes, title){
    data<-results[1:(min(nrow(data),10)),]
    geneRatio<-data$nbhits/length(listGenes)
    pal<-colorRampPalette(c("red", "blue"))
    col<-pal(length(unique(data$p.adj)))[sapply(data$p.adj, function(p){which(unique(data$p.adj==p))})]
    par(mar=c(5,max(sapply(rownames(data), nchar))/1.5, 4,6), xpd=TRUE)
    plot(geneRatio, 1:nrow(data)+1, col=col, pch=19, yaxt="n", ylab="", xlab="Gene Ratio", cex=log10(data$nbhits)+1, main=title)
    axis(2, at=1:nrow(data)+1, label=rownames(data), las=2)
    legend("topright", inset=c(-0.14, -0.04), 
           legend =unique(c(min(data$p.adj),median(data$p.adj), max(data$p.adj))), 
           pch = 19, 
           col = pal(length(unique(data$p.adj)))[sapply(unique(c(min(data$p.adj),median(data$p.adj), max(data$p.adj))), function(level){
             result<-which(unique(data$p.adj)==level)
             if(length(result)==0)return(0)
             else return(result[1])
           })],
    title="P value", bty="n")
    legend("topright", inset=c(-0.14, 0.5), legend=unique(data$nbhits), pch=19, col="gray", pt.cex=log10(unique(data$nbhits))+1, bty="n", title="Hit\n number")
  }
  
  # Select geneSet -------------------------------------------------------------------------------------------------------------------------------------------
  geneSetDomain<-reactive({
    return(geneSubset("SEA", data(), input$geneIdDomain, input$criterionDomain, input$setSelectDomain, input$upperSelectDomain))
  })
  
  # Run test -------------------------------------------------------------------------------------------------------------------------------------------------
  resultDomain<-eventReactive(input$runDomain,{

    # Annotate samples
    samples <-getBM(attributes = c('ensembl_gene_id', 'family'),
                    filters = 'ensembl_gene_id', #implémenter conversion
                    values = geneSetDomain(), 
                    mart = mart)
    samples<-samples[,2]
    # Get all possible values
    universe<-getBM(attributes = c('ensembl_gene_id', 'family'),
                    mart = mart)
    universe<-universe[,2]
    
    # Run test
    return(enrich.domains(samples, unique(samples), universe, p.adj.method=input$correctionDomain))
    
  })
  
  # Dotplot ------------------------------------------------------------------------------------------------------------------
  output$dotplotDomain<-renderPlot({
    if(! is.null(resultDomain())) {
      dotplot.domains(resultDomain(), geneSetDomain(), title="Top 10 enriched domains \n (Panther)")
    }
  })
  
})