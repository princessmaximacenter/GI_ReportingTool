# Purpose of the script: Collecting information of a specific candidate 
# genetic interaction
#
# Author(s): Denise Kersjes
# Date of creation:  28  August 2020
# Date of last edit: 30  November 2020

script.version <- "1.0"


### Installation
source("setup.R")
source("functions.R")


### Start the interactive application

# Define UI for application that draws a histogram
ui <- fluidPage(
    
    # Create a first panel to for setting the title
    wellPanel(fluidRow(
      column(
        width = 11,
        
        # Set the title at the top in the middle
        align = "center",
        titlePanel("Candidate genetic interactions in childhood cancer")
      )
    )),
    
    # Create a second panel to select a candidate
    fluidRow(
      
      # Column at the left side of the browser
      column(
        width = 4,
        
        # Give some information about selecting the candidate
        p(tags$br(),
          "Select a cancer type and a candidate of interest to get more \
          information about the candidate genetic interaction. Please, note \
          the gene names are in aplhabetic order. Co-occuring candidates are \
          indicated as 'CO' and mutually exclusive as 'ME'.",
          style="text-align: justify;"),
        
        # Insert a drop down menu to select the cancer type
        selectInput(inputId = "cancertype",
                    label = "Select the cancer type",
                    choices = sort(unique(cand.df$cancer.type))),

        # Insert a drop down menu to select the candidate
        selectInput(inputId = "candidate",
                    label = "Select the candidate",
                    choices = NULL)
      ),
      
      # Column at the right side of the browser
      column(
        width = 8,
        align = "center",
        
        # Insert the overview map of the candidates
        p(tags$br()),
        img(src="GI_map.png",
            width = "90%"),
        p(tags$br())
      )
    ),
    
    # Create a third panel for the candidate information
    wellPanel(fluidRow(
      column(
        width = 11,
        align = "justify",
        
        # Set the tab panels
        tabsetPanel(
          tabPanel("Variant summary of cancer type", htmlOutput("infoSummary"), 
                   plotOutput("summary")),
          tabPanel("Lolliplot of protein changes", htmlOutput("infoLolliplot"),
                   plotOutput("lolliplot1"), plotOutput("lolliplot2")), 
          tabPanel("Heatmap of genomic alterations", htmlOutput("infoOncoprint"), 
                   plotOutput('oncoprint')),
          tabPanel("Variant effect predictor", htmlOutput("infoVep"), 
                   dataTableOutput('vep')))
      ))
    ),
    
    # Add the version of the application
    fluidRow(
      
      # Column at the left side of the browser
      column(
        width = 12,
        align = "right",
        p(paste0("version ", script.version))
      ))
    )
  
# Define server logic required to get analysis results
server <- function(input, output, session) {
    
    # Set the correct candidates for each cancer type
    observe({
      # Get the selected cancer type 
      ct <- input$cancertype
      
      # Specify the candidates per cancer type
      ct.cand <- dplyr::filter(cand.df, cand.df$cancer.type == ct)
      ct.cand.visible <- base::paste0(ct.cand$gene1, " - ", 
                                      ct.cand$gene2, " (", 
                                      base::toupper(ct.cand$co.me), ", ", 
                                      base::toupper(ct.cand$dataset) , ")")
        
      # Update candidates when changing cancer type
      updateSelectInput(session, "candidate",
                        choices = ct.cand.visible
      )
    })
    
    # Get the gene names and the dataset of the selected candidate
    selected_candidate <- reactive({
      
      splitted <- BiocGenerics::unlist(base::strsplit(
        x = input$candidate, split = " "))
      gene1 <- splitted[1]
      gene2 <- splitted[3]
      dataset <- stringr::str_sub(splitted[5], end = -2)
    
      return(list(gene1 = gene1, gene2 = gene2, dataset = dataset))
    })
  
    # Get the stable gene names of the selected candidate
    ensemble_id <- reactive({
      
      # Get the Ensembl stable ID of gene 1
      gene1.df <- dplyr::filter(vep.df, 
                                vep.df$SYMBOL == selected_candidate()$gene1)
      gene1.ens.id <- IRanges::unique(gene1.df$Gene)
 
      # Get the Ensembl stable ID of gene 2
      gene2.df <- dplyr::filter(vep.df, 
                                vep.df$SYMBOL == selected_candidate()$gene2)
      gene2.ens.id <- IRanges::unique(gene2.df$Gene)
      
      return(list(stable.gene1 = gene1.ens.id, stable.gene2 = gene2.ens.id))
    })
    
    # Merge the mutation information file with the VEP output file
    create_maf <- reactive({
      
      # Select the cancer type 
      if (input$cancertype != "PAN") {
        maf.df <- maf.df[IRanges::which(maf.df$ct == input$cancertype), ] 
      }
      # Select the data set 
      maf.df <- maf.df[which(
        maf.df$dataset == base::tolower(selected_candidate()$dataset)), ] 
      
      # Read the final data frame as a MAF object
      maf.object <- read.maf(maf = maf.df)
      
      return(maf.object)
    })
    
    # Create a variant summary of the cancer type
    create_summary <- reactive({
      plotmafSummary(maf = create_maf(), 
                     rmOutlier = FALSE, 
                     dashboard = TRUE)
    })
    
    # Create lolliplots of amino acid changes of the candidate
    create_lolliplot_1 <- reactive({
      lollipopPlot(maf = create_maf(),
                   gene = c(selected_candidate()$gene1),
                   AACol = "HGVSp_short",
                   showDomainLabel = FALSE,
                   labelPos = "all")
    })
    create_lolliplot_2 <- reactive({
      lollipopPlot(maf = create_maf(),
                   gene = c(selected_candidate()$gene2),
                   AACol = "HGVSp_short",
                   showDomainLabel = FALSE,
                   labelPos = "all")
    })
    
    # Create a heatmap of genomic alterations of the candidate
    create_oncoprint <- reactive({
      oncoplot(maf = create_maf(), 
               genes = c(as.character(selected_candidate()$gene1), 
                         as.character(selected_candidate()$gene2)))
    })
    
    # Extract information of Ensembl Variant Effect Predictor (VEP)
    vep_info <- reactive({
      
      # Extract transcipts of the candidate genes 
      vep.info.genes <- merged.df[merged.df$Gene %in% ensemble_id(), ]
      
      # Order the list based on gene name
      transcripts.genes <- vep.info.genes[GenomicRanges::order(
        vep.info.genes$Gene),] 
      
      # Extract a desired subset of the data
      columns.of.interest <-
        c("Hugo_Symbol", "Allele", "Variant_Type", "Consequence", 
          "BIOTYPE", "HGVSp_short", "IMPACT", "SIFT", "PolyPhen", "Condel", 
          "Existing_variation")
      sub.info.transcripts <- transcripts.genes[columns.of.interest]
      
      return(sub.info.transcripts)
    })
    
    # Set information about the variant summary
    output$infoSummary <- renderText(
      paste("<p><br> A summary of the pediatric", input$cancertype, "of the", 
            selected_candidate()$dataset,
            "dataset is given below. The distribution of variant classification, \
            variant type, and the mutational change is visualized at the top. \
            The boxplot at the bottom gives a more detailed distribution of the \
            variant classification. The stacked barplot shows the number of \
            variants in each sample and the variant distribution of the top \
            mutated genes in", input$cancertype, "is shown in the last plot. \
            <p><br>", sep = " ")
      )
    # Display the variant summary of the cancer type
    output$summary <- renderPlot(
      create_summary()
    )
    
    # Set information about the lolliplot
    output$infoLolliplot <- renderText(
      paste("<p><br> A lolliplot visualize which protein changes occur in the \
            data and where they are located at the gene. Lolliplots of the \
            candidate genes ", selected_candidate()$gene1, "and ", 
            selected_candidate()$gene2, "are shown below. <p><br>", sep = " ")
    )
    # Display the lolliplots of amino acid changes of the candidate
    output$lolliplot1 <- renderPlot(
      create_lolliplot_1()
    )
    output$lolliplot2 <- renderPlot(
      create_lolliplot_2()
    )
    
    # Set information about the oncoprint plot
    output$infoOncoprint <- renderText(
      paste("<p><br> OncoPrint is a way to visualize genomic alterations by \
            heatmap. Every column is a sample in the data. A colored bar \
            indicate a mutation in the on the right mentioned gene. The total \
            number of mutations in a sample or gene are visualized in the top \
            barplot and side barplot respectively. <p><br>", sep = " ")
      )
    # Display the heatmap of genomic alterations of the candidate
    output$oncoprint <- renderPlot(
      create_oncoprint()
    )
    
    # Set information about the VEP table
    output$infoVep <- renderText(
      paste("<p><br> Mutations in", selected_candidate()$gene1, "and", 
            selected_candidate()$gene2, "are analysed in Ensembl Variant Effect \
            Predictor (VEP) GRCh37. The effect of the transcript of each sample \
            is shown in the table below. <p><br>", sep = " ")
      )
    # Display the VEP information table
    output$vep <- renderDataTable(
      vep_info(),
      options = list(pageLength = 10)
    )
  }
  
# Run the application 
shinyApp(ui = ui, server = server)
