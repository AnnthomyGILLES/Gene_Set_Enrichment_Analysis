library(shiny)
library(shinydashboard)
library(clusterProfiler)
library(AnnotationDbi)
library(GO.db)
library(topGO)
library(ReactomePA)
library(biomaRt) 
library(DT)
library(shinythemes)
library(shinyjs)
library(shinycssloaders)
library(DOSE)
library(shinyBS)

# Define UI for application that draws a histogram
shinyUI(
  tagList(useShinyjs(),
  navbarPage("Make Shiny Great Again! v1.0",id="enrichAnno",theme = shinytheme("flatly"),
             tabPanel('DATA', verbatimTextOutput('fileInfo'),
                      pageWithSidebar(
                        headerPanel("Import your Dataset"),
                        # sidebarLayout(
                        sidebarPanel(
                          fileInput('inFile', 'Choose CSV File', accept = c('text/csv', 'text/comma-separated-values,text/plain', '.csv')),
                          checkboxInput('header', 'Header', TRUE),
                          radioButtons('sep', 'Separator', c(Comma = ',', Semicolon = ';', Tab = '\t'), '\t'),
                          radioButtons('quote', 'Quote', c(None = '', 'Double Quote'='"', 'Single Quote' = "'"),'"'),
                          tags$hr(), 
                          width = 3
                        ),
                        mainPanel(
                          DT::dataTableOutput("eventTable") 
                        )
                      )),
             tabPanel("GO",value = "GO",
                      # sidebarLayout(
                      pageWithSidebar(
                        
                        # Application title
                        headerPanel("Gene Ontology"),
                        
                        sidebarPanel(

                          # # Select columns and gene set
                          tipify(uiOutput(outputId="IdSelector"),"Name of column containing gene IDs",placement = "right"),
                          tipify(uiOutput(outputId="OrderSelector"),"Name of column containing criterion for gene ordering/selection",placement = "right"),
                          tipify(uiOutput(outputId="upper"),"Select genes with criterion higher/lower than cutoff, default (unticked) is lower",placement = "right"),
                          tipify(uiOutput(outputId="setCutoff"),"Criterion cutoff for gene selection",placement = "right"),
                          
                          # Annotation options
                          tipify(uiOutput(outputId = "annOption2"),"filter GO enriched result at specific level",placement = "right"),
                          uiOutput(outputId = "annOption3"),
                          # 
                          # # Test selection
                          tipify(selectInput(inputId = "test", label = "Select test to run :",
                                      choices = list( "SEA", "GSEA")),"SEA : Subset Enrichment Analysis, uses Fisher exact test to estimate probability of annotations to occur in this subset. GSEA : Gene Set Enrichment Analysis uses random permutation test to estimate probability of annotations to be associated with higher/lower criterion values",placement = "right"),
                          
                          # # Test options
                          tipify(uiOutput(outputId = "testpCutoff"),"Maximum p value of conserved results",placement = "right"),
                          tipify(uiOutput(outputId = "testqCutoff"),"Maximum q value of conserved results",placement = "right"),
                          tipify(uiOutput(outputId = "testSetSize"),"Minimum and maximum number of genes for a set to be tested",placement = "right"),
                          tipify(uiOutput(outputId = "testCorrection"),"bonferonni= Bonferonni correction (more stringent), BH = Benjamini Hochberg s procedure (less stringent), none = uncorrected",placement = "right"),
                          tipify(checkboxInput("donum1", "BP", value = T),"Biological Process",placement = "right"),
                          tipify(checkboxInput("donum2", "CC", value = T),"Cellular Component",placement = "right"),
                          tipify(checkboxInput("donum3", "MF", value = T),"Molecular Function",placement = "right"),
                          actionButton(inputId="plotWholeGO", label="Barplot"),
                          tags$div(
                            tags$br()
                          ),
                          actionButton(inputId="runGO", label="Run test"),
                          tags$div(
                            tags$br()
                          ),
                          # actionButton(inputId="goGraph",label="GO Graph"),
                          width = 3
                        ),
                        mainPanel(
                          tabsetPanel(
                            tabPanel("BARPLOT", value = "BARPLOT",
                                     bsCollapse(id = "collapseExample_gsea",open = "GROUPGO_SEA_BP",
                                                bsCollapsePanel("Biological Process",conditionalPanel("input.donum1==true && input.plotWholeGO>0",withSpinner(plotOutput(outputId = "wholeGOSetBP"))),value = "GROUPGO_SEA_BP", style="primary"),
                                                bsCollapsePanel("Cellular Component",conditionalPanel("input.donum2==true && input.plotWholeGO>0",withSpinner(plotOutput(outputId = "wholeGOSetCC"))),value = "GROUPGO_SEA_CC", style="info"),
                                                bsCollapsePanel("Molecular Function",conditionalPanel("input.donum3==true && input.plotWholeGO>0",withSpinner(plotOutput(outputId = "wholeGOSetMF"))),value = "GROUPGO_SEA_MF", style="default"),multiple=TRUE
                                     )
                            ),
                            tabPanel("DOTPLOT", value = "DOTPLOT",
                                     bsCollapse(id = "collapseExample",open = "dotplot_BP",
                                                bsCollapsePanel("Dotplot BP",conditionalPanel("input.donum1==true  && input.runGO>0",withSpinner(plotOutput(outputId = "dotplotBP"))),value = "dotplot_BP", style="primary"),
                                                bsCollapsePanel("Dotplot CC",conditionalPanel("input.donum2==true  && input.runGO>0",withSpinner(plotOutput(outputId = "dotplotCC"))),value = "dotplot_CC", style="info"),
                                                bsCollapsePanel("Dotplot MF",conditionalPanel("input.donum3==true  && input.runGO>0",withSpinner(plotOutput(outputId = "dotplotMF"))),value = "dotplot_MF", style="default"),multiple = TRUE
                                     )
                            ),
                            tabPanel("GOGRAPH", value = "GOGRAPH",
                                     bsCollapse(id = "collapseExample",open = "plotGOgraphBP",
                                                bsCollapsePanel("GO graph BP",conditionalPanel("input.donum1==true  && input.runGO>0 && input.test=='SEA'",withSpinner(plotOutput(outputId = "plotGOgraphBP"))),value = "plotGOgraphBP", style="primary"),
                                                bsCollapsePanel("GO graph CC",conditionalPanel("input.donum2==true  && input.runGO>0 && input.test=='SEA'",withSpinner(plotOutput(outputId = "plotGOgraphCC"))),value = "plotGOgraphCC", style="info"),
                                                bsCollapsePanel("GO graph MF",conditionalPanel("input.donum3==true  && input.runGO>0 && input.test=='SEA'",withSpinner(plotOutput(outputId = "plotGOgraphMF"))),value = "plotGOgraphMF", style="default"),multiple = TRUE
                                                
                                     )
                            )
                          )
                        )
   
                      )
                  ),
             tabPanel("PATHWAY",value = "Pathway",
                      # Application title
                      headerPanel("Pathways"),
                      sidebarLayout(
                        sidebarPanel(
                          # Select columns and gene set
                          #Id selection
                          tipify(uiOutput(outputId="IdSelectorPathway"),"Name of column containing gene IDs",placement = "right"),
                          # Ordering criterion selection
                          tipify(uiOutput(outputId="OrderSelectorPathway"),"Name of column containing criterion for gene ordering/selection",placement = "right"),
                          # Select gene set by upper values
                          tipify(uiOutput(outputId="upperPathway"),"Select genes with criterion higher/lower than cutoff, default (unticked) is lower",placement = "right"),
                          # Select selection threshold
                          tipify(uiOutput(outputId="setCutoffPathway"),"Criterion cutoff for gene selection",placement = "right"),

                          # Annotation options
                          tipify(uiOutput(outputId="annOption1Pathway"),"Reference database from which annotations are extracted",placement = "right"),
                          
                          # Test selection
                          tipify(selectInput(inputId = "testPathway", label = "Select test to run :",
                                      choices = list( "SEA", "GSEA")),"SEA : Subset Enrichment Analysis, uses Fisher's exact test to estimate probability of annotations to occur in this subset \n GSEA : Gene Set Enrichment Analysis uses random permutation test to estimate probability of annotations to be associated with higher/lower criterion values",placement = "right"),

                          # Test options
                          tipify(uiOutput(outputId = "testpCutoffPathway"),"Maximum p value of conserved results",placement = "right"),
                          tipify(uiOutput(outputId = "testqCutoffPathway"),"Maximum q value of conserved results",placement = "right"),
                          tipify(uiOutput(outputId = "testSetSizePathway"),"Minimum and maximum number of genes for a set to be tested",placement = "right"),
                          tipify(uiOutput(outputId = "testCorrectionPathway"),"bonferonni= Bonferonni correction (more stringent), BH = Benjamini Hochberg s procedure (less stringent), none = uncorrected",placement = "right"),
                          

                          # Display result
                          uiOutput(outputId = "pathwayResultIds"),
                          
                          #Run test
                          uiOutput(outputId = "runTestPathway")
                          # plotOutput(outputId = "enrichMapPathway"),
                          # plotOutput(outputId = "pathwayPlotPathway")
                        ),
                        mainPanel(
                          tabsetPanel(
                            tabPanel("RESULTS", value = "Pathway",
                                     conditionalPanel("input.testPathway=='SEA'",
                                     bsCollapse(id = "collapsePathwaySEA",
                                                bsCollapsePanel("Barplot",conditionalPanel("input.runPathway>0",withSpinner(plotOutput(outputId = "barplotPathway"))),value = "barplotPathway", style="primary"),
                                                bsCollapsePanel("Dotplot",conditionalPanel("input.runPathway>0",withSpinner(plotOutput(outputId = "dotplotPathway"))),value = "dotplotPathway", style="info"),
                                                bsCollapsePanel("EnrichMap",conditionalPanel("input.runPathway>0",withSpinner(plotOutput(outputId = "enrichMapPathway"))),value = "enrichMapPathway", style="success"),
                                                bsCollapsePanel("Pathview",conditionalPanel("input.runPathway>0 && input.database=='KEGG' ",withSpinner(plotOutput(outputId = "pathwayPlot"))),value = "view Pathway", style="default"),multiple = TRUE
                                                # bsCollapsePanel("enrichMap",conditionalPanel("input$runGO",withSpinner(plotOutput(outputId = "GOSetMF"))),value = "GROUPGO_GSEA_MF", style="warning")
                                     )),
                                     conditionalPanel("input.testPathway=='GSEA'",
                                                      bsCollapse(id = "collapsePathwayGSEA",open = c("Dotplot","EnrichMap"),
                                                                 bsCollapsePanel("Dotplot",conditionalPanel("input.runPathway>0",withSpinner(plotOutput(outputId = "dotplotPathwayGSEA"))),value = "dotplotPathwayGSEA", style="info"),
                                                                 bsCollapsePanel("EnrichMap",conditionalPanel("input.runPathway>0",withSpinner(plotOutput(outputId = "enrichMapPathwayGSEA"))),value = "enrichMapPathwayGSEA", style="success"),
                                                                 bsCollapsePanel("Pathview",conditionalPanel("input.runPathway>0 && input.database=='KEGG'",withSpinner(plotOutput(outputId = "pathwayPlotGSEA"))),value = "view Pathway", style="default"),multiple = TRUE
                                                                 
                                                      ))
                            )
                          )
                        )
                      )
             ),
             tabPanel("DOMAIN",value = "Domain",
                      # Application title
                      headerPanel("Domains"),
                      sidebarLayout(
                        sidebarPanel(
                          
                          # Select columns and gene set
                          #Id selection
                          tipify(uiOutput(outputId="IdSelectorDomain"),"Name of column containing gene IDs",placement = "right"),
                          # Ordering criterion selection
                          tipify(uiOutput(outputId="OrderSelectorDomain"),"Name of column containing criterion for gene ordering/selection",placement = "right"),
                          # Select gene set by upper values
                          tipify(uiOutput(outputId="upperDomain"),"Select genes with criterion higher/lower than cutoff, default (unticked) is lower",placement = "right"),
                          # Select selection threshold
                          tipify(uiOutput(outputId="setCutoffDomain"),"Criterion cutoff for gene selection",placement = "right"),

                          #tableOutput(outputId="contents"),

                          # Test options
                          tipify(uiOutput(outputId = "testCorrectionDomain"),"bonferonni= Bonferonni's correction (more stringent) \n BH = Benjamini Hochberg's procedure (less stringent)",placement = "right"),
                          
                          # Display result

                          #Run test
                          uiOutput(outputId = "runTestDomain")
                          # plotOutput(outputId = "dotplot")
                          ),
                        mainPanel(
                          bsCollapse(id = "collapseExampleDomain",open = "dotplotDOmain",
                                     bsCollapsePanel("Dotplot BP",conditionalPanel("input$runTestDomain",withSpinner(plotOutput(outputId = "dotplotDomain"))),value = "dotplotDOmain", style="success")
                          )
                        )
                      )
             ),
             tabPanel("ABOUT",value = "About",
               div(
                 h1("License"),
                 p("This software was developed by GILLES Annthomy and BRELURUT Geoffray and is available freely."),
                p("The source code is available on GitHub.For any comments or questions, please contact annthomy.gilles@etu.univ-rouen.fr or geoffray.brelurut@etu.univ-rouen.fr")
               )
              ),
              collapsible = T
            )
  )
)