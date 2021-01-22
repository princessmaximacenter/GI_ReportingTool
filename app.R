# Purpose of the script: Collecting information of a specific candidate 
# genetic interaction
#
# Author(s): Denise Kersjes
# Date of creation:  28  August 2020
# Date of last edit: 22  January 2021

script.version <- "1.1"


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
          information about the candidate genetic interaction. Note that \
          the gene names are in aplhabetic order. Co-occuring candidates are \
          indicated as 'CO' and mutually exclusive as 'ME'.",
          style="text-align: justify;"),
        
        # Insert a drop down menu to select the cancer type
        selectInput(inputId = "cancertype",
                    label = "Select the cancer type",
                    choices = sort(unique(cand.df$drop_list))),

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
                   radioButtons(inputId = "come", label = "Variants to show",
                                choices = list("all variants", "only CO variants", 
                                           "only ME variants")),
                   textOutput("radio"),
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
      ct.cand <- dplyr::filter(cand.df, cand.df$drop_list == ct)
      ct.cand.visible <- base::paste0(ct.cand$gene1, " - ", 
                                      ct.cand$gene2, " (", 
                                      base::toupper(ct.cand$co.me), ")")
        
      # Update candidates when changing cancer type
      updateSelectInput(session, "candidate",
                        choices = ct.cand.visible
      )
    })
  

    # Get the gene names and the dataset of the selected candidate
    selected_candidate <- reactive({
      
      splitted_drop1 <- BiocGenerics::unlist(base::strsplit(
        x = input$cancertype, split = " "))
      cancertype <- splitted_drop1[1]
      dataset <- stringr::str_sub(splitted_drop1[2], start = 2, end = -2)
      
      splitted_drop2 <- BiocGenerics::unlist(base::strsplit(
        x = input$candidate, split = " "))
      gene1 <- splitted_drop2[1]
      gene2 <- splitted_drop2[3]

      return(list(gene1 = gene1, gene2 = gene2, dataset = dataset, ct = cancertype))
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
    
    # Get patients that have a mutation in both genes of the candidate
    co_patients <- reactive({
      
      # Get subpart of data frame wiht only the candidate genes
      merged.df.genes <- merged.df[merged.df$Gene %in% ensemble_id(), ]
      
      # Get the patients of the first gene
      gene1.df <- merged.df.genes[which(merged.df.genes["Gene"] == 
                                         ensemble_id()$stable.gene1), ]
      patients.gene1 <- unique(gene1.df$Tumor_Sample_Barcode)

      # Get the patients of the second gene
      gene2.df <- merged.df.genes[which(merged.df.genes["Gene"] == 
                                         ensemble_id()$stable.gene2), ]
      patients.gene2 <- unique(gene2.df$Tumor_Sample_Barcode)

      # Get the overlapping patients
      overlap_patients <- patients.gene1[patients.gene1 %in% patients.gene2]
      
      return(list(ids = overlap_patients))
    })
    
    # Merge the mutation information file with the VEP output file
    create_maf <- reactive({
      
      # Select the cancer type 
      if (selected_candidate()$ct != "PAN") {
        maf.df <- maf.df[IRanges::which(maf.df$ct == selected_candidate()$ct), ] 
      }
      # Select the data set 
      maf.df <- maf.df[which(
        maf.df$dataset == base::tolower(selected_candidate()$dataset)), ] 
      
      # Indicate if the patients of the candidates have a mutation in both genes
      maf.df <- mutate(maf.df, CO_ME = ifelse(
        Tumor_Sample_Barcode %in% co_patients()$ids, "CO", "ME"))
      
      # Create subset if only CO or ME candidates should be shown
      if (input$come == "only CO variants") {
        maf.df <- maf.df[IRanges::which(maf.df$Tumor_Sample_Barcode 
                                        %in% co_patients()$ids), ]
      } else if (input$come == "only ME variants") {
        maf.df <- maf.df[IRanges::which(!maf.df$Tumor_Sample_Barcode 
                                        %in% co_patients()$ids), ]
      }
        
      # Read the final data frame as a MAF object
      maf.object <- read.maf(maf = maf.df)
      
      return(maf.object)
    })
    
    # Create a variant summary of the cancer type
    create_summary <- reactive({
      # Update first the radio button to include alle variants
      updateRadioButtons(session, inputId = "come", 
                         choices = list("all variants", "only CO variants", 
                                        "only ME variants"))
      plotmafSummary(maf = create_maf(), 
                     rmOutlier = FALSE, 
                     dashboard = TRUE,
                     titleSize = 2,
                     fs = 1.4)
    })
    
    # Create lolliplots of amino acid changes of the candidate
    create_lolliplot_1 <- reactive({
      lollipopPlot(maf = create_maf(),
                   gene = c(selected_candidate()$gene1),
                   AACol = "HGVSp_short",
                   showDomainLabel = FALSE,
                   labelPos = "all",
                   labPosSize = 4,
                   legendTxtSize = 15,
                   axisTextSize = c(15, 15),
                   titleSize = c(20, 18),
                   pointSize = 2)
    })
    create_lolliplot_2 <- reactive({
      lollipopPlot(maf = create_maf(),
                   gene = c(selected_candidate()$gene2),
                   AACol = "HGVSp_short",
                   showDomainLabel = FALSE,
                   labelPos = "all",
                   labPosSize = 4,
                   legendTxtSize = 15,
                   axisTextSize = c(15, 15),
                   titleSize = c(20, 18),
                   pointSize = 2)
    })
    
    # Create a heatmap of genomic alterations of the candidate
    create_oncoprint <- reactive({
      # Update first the radio button to include alle variants
      updateRadioButtons(session, inputId = "come", 
                         choices = list("all variants", "only CO variants", 
                                                      "only ME variants"))
      oncoplot(maf = create_maf(), 
               genes = c(as.character(selected_candidate()$gene1), 
                         as.character(selected_candidate()$gene2)))
    })
    
    # Extract information of Ensembl Variant Effect Predictor (VEP)
    vep_info <- reactive({
      
      # Extract transcipts of the candidate genes 
      vep.info.genes <- merged.df[merged.df$Gene %in% ensemble_id(), ]
      
      # Extract only genes belonging to the selected cancer type and dataset
      vep.info.genes <- vep.info.genes[(vep.info.genes$ct == 
                          selected_candidate()$ct & vep.info.genes$dataset == 
                            tolower(selected_candidate()$dataset)), ]

      # Check which patients have a mutation in both genes
      vep.info.genes <- mutate(vep.info.genes, CO_ME = ifelse(
        Tumor_Sample_Barcode %in% co_patients()$ids, "CO", "ME"))
      
      # Order the list based on gene name
      transcripts.genes <- vep.info.genes[GenomicRanges::order(
        vep.info.genes$Gene),] 
      
      # Extract a desired subset of the data
      columns.of.interest <-
        c("Hugo_Symbol", "CO_ME", "Chromosome", "Start_Position", "End_Position", 
          "Allele", "Variant_Type", "Consequence", "BIOTYPE", "HGVSp_short", 
          "IMPACT", "SIFT", "PolyPhen", "Condel", "Existing_variation")
      sub.info.transcripts <- transcripts.genes[columns.of.interest]
      
      return(sub.info.transcripts)
    })
    
    # Set information about the variant summary
    output$infoSummary <- renderText(
      paste("<p><br> A summary of pediatric <b>", selected_candidate()$ct, 
            "</b> of the <b>", selected_candidate()$dataset, "</b> dataset \
            is given below. The distribution of variant classification, \
            variant type, and the mutational change is visualized at the top. \
            The stacked barplot at the bottom shows the number of variants in \
            each sample. The boxplot gives a more detailed distribution of the \
            variant classification. The variant distribution of the top \
            mutated genes in", selected_candidate()$ct, "is shown in the last \
            plot. <p><br>", sep = " ")
      )
    # Display the variant summary of the cancer type
    output$summary <- renderPlot(
      create_summary()
    )
    
    # Set information about the lolliplot
    output$infoLolliplot <- renderText(
      paste("<p><br> A lolliplot visualizes which protein changes occur in the \
            data and where they are located at the gene. Lolliplots of the \
            candidate genes <b>", selected_candidate()$gene1, "</b> and <b>", 
            selected_candidate()$gene2, "</b> are shown below. If the error \
            of <em>'No non-synonymous mutations found'</em> pops up after \
            selecting only CO or ME variants, this means the candidate does \
            not have either CO or ME variants. <p><br>", sep = " ")
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
      paste("<p><br> An oncoPrint is a way to visualize genomic alterations \
            with a heatmap. Every column is a sample in the data. A colored bar \
            indicates a mutation in the gene mentioned on the right. The total \
            number of mutations in a sample or gene is visualized in the top \
            barplot and side barplot respectively. <p><br>", sep = " ")
      )
    # Display the heatmap of genomic alterations of the candidate
    output$oncoprint <- renderPlot(
      create_oncoprint()
    )
    
    # Set information about the VEP table
    output$infoVep <- renderText(
      paste("<p><br> Mutations in <b>", selected_candidate()$gene1, "</b> and \
            <b>", selected_candidate()$gene2, "</b> are analysed in Ensembl \
            Variant Effect Predictor (VEP) GRCh37. The effect of the transcript \
            of each sample is shown in the table below. <p><br>", sep = " ")
      )
    # Display the VEP information table
    output$vep <- renderDataTable(
      vep_info(),
      options = list(pageLength = 10)
    )
  }
  
# Run the application 
shinyApp(ui = ui, server = server)
