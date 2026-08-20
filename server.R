library(shiny)
library(Seurat)
library(ggplot2)
library(dplyr)

human <- readRDS("data/Human_DietSeurat.rds")
mouse <- readRDS("data/Mouse_DietSeurat.rds")
human_to_mouse <- readRDS("data/human_to_mouse.rds")
mouse_to_human <- readRDS("data/mouse_to_human.rds")

human_color <- "#B3D9FF"
mouse_color <- "#FFCCCC"

get_gene_pair <- function(gene) {
  
  human_gene <- NA_character_
  mouse_gene <- NA_character_
  
  if (gene %in% rownames(human)) {
    human_gene <- gene
    mouse_gene <- unname(human_to_mouse[gene])
  }
  
  if (gene %in% rownames(mouse)) {
    mouse_gene <- gene
    human_gene <- unname(mouse_to_human[gene])
  }
  
  list(
    human_gene = human_gene,
    mouse_gene = mouse_gene
  )
}


make_vln <- function(obj, gene, fill_color, species_label) {
  
  group_col <- "clusters"
  
  validate(
    need(!is.na(gene), paste0("Undetectable expression, or no ", species_label, " homolog found."))
  )
  
  df <- FetchData(
    obj,
    vars = c(gene, group_col)
  )
  
  colnames(df) <- c("expression", "group")
  
  ggplot(df, aes(x = group, y = expression)) +
    geom_violin(
      fill = fill_color,
      color = "black",
      trim = TRUE,
      scale = "width"
    ) +
    geom_boxplot(
      width = 0.05,
      fill = "white",
      color = "black",
      outlier.size = 0.2
    ) +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 18),
      axis.text.x = element_text(size = 16),
      axis.text.y = element_text(size = 16)
    ) +
    labs(
      title = gene,
      x = NULL,
      y = "Lognorm Expression"
    )
}


server <- function(input, output, session) {
  
  human_genes <- setNames(rownames(human), paste0(rownames(human), " - human"))
  mouse_genes <- setNames(rownames(mouse), paste0(rownames(mouse), " - mouse"))
  
  gene_choices <- sort(c(
    human_genes,
    mouse_genes
  ))

  priority_genes <- gene_choices[nchar(gene_choices) <= 3]
  
  gene_choices <- c(
    priority_genes,
    gene_choices[!gene_choices %in% priority_genes]
  )
  
  updateSelectizeInput(
    session = session,
    inputId = "gene",
    choices = c("", gene_choices),
    selected = "",
    server = TRUE,
    options = list(
      placeholder = "GENE SYMBOL",
      dropdownParent = "body",
      openOnFocus = FALSE,
      highlight = FALSE,
      maxOptions = 100,
      minLength = 2,
      
      ## matching "beginWith"
      score = I("function(search) {
        return function(item) {
          // Strict server-side alignment check using '^' anchoring
          var regex = new RegExp('^' + search, 'i');
          return item.label.match(regex) ? 1 : 0;
        };
      }")
      
    )
  )

  gene_pair <- reactive({
    req(input$gene)
    get_gene_pair(input$gene)
  })
  
  output$human_vln <- renderPlot({
    req(!is.null(input$gene),
        nzchar(input$gene),
        cancelOutput = TRUE
    )
    
    gp <- gene_pair()
    
    make_vln(
      obj = human,
      gene = gp$human_gene,
      fill_color = human_color,
      species_label = "human"
    )
  })
  
  output$mouse_vln <- renderPlot({
    req(!is.null(input$gene),
        nzchar(input$gene),
        cancelOutput = TRUE
    )
    
    gp <- gene_pair()
    
    make_vln(
      obj = mouse,
      gene = gp$mouse_gene,
      fill_color = mouse_color,
      species_label = "mouse"
    )
  })
  
  selected_source <- reactiveVal("source_Grinberg2026")
  
  observeEvent(input$source_Grinberg2026, {
    selected_source("source_Grinberg2026")
  })
  
  # observeEvent(input$source_mouse_microglia, {
  #   selected_source("source_mouse_microglia")
  # })
  # 
  # observeEvent(input$source_mouse_astrocyte, {
  #   selected_source("source_mouse_astrocyte")
  # })
  
}