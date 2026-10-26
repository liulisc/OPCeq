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


get_cluster_levels <- function(obj) {
  md <- tryCatch(slot(obj, "meta.data"), error = function(e) NULL)
  if (is.null(md) || !"clusters" %in% colnames(md)) {
    return(character(0))
  }
  cl <- md$clusters
  if (is.factor(cl)) {
    levels(cl)
  } else {
    sort(unique(as.character(cl)))
  }
}

fetch_gene_df <- function(obj, gene) {
  df <- FetchData(obj, vars = c(gene, "clusters"))
  colnames(df) <- c("expression", "group")
  df$group <- factor(as.character(df$group), levels = get_cluster_levels(obj))
  df
}

sig_status <- function(p_value) {
  if (is.na(p_value) || p_value >= 0.05) {
    "Insignificant"
  } else {
    "Significant"
  }
}

p_to_label <- function(p_value) {
  if (is.na(p_value) || p_value >= 0.05) {
    return("N/S")
  }
  if (p_value < 0.001) {
    return("***")
  }
  if (p_value < 0.01) {
    return("**")
  }
  "*"
}

safe_t_test <- function(a, b) {
  tryCatch(
    t.test(a, b),
    error = function(e) NULL
  )
}

safe_t_p <- function(a, b) {
  tt <- safe_t_test(a, b)
  if (is.null(tt)) {
    NA_real_
  } else {
    unname(tt$p.value)
  }
}

add_ttest_star <- function(p, df, g1, g2, y_limit = NULL) {
  x_levels <- levels(df$group)
  x1 <- match(g1, x_levels)
  x2 <- match(g2, x_levels)
  if (is.na(x1) || is.na(x2)) {
    return(p)
  }
  
  a <- df$expression[as.character(df$group) == g1]
  b <- df$expression[as.character(df$group) == g2]
  if (length(a) < 2 || length(b) < 2) {
    return(p)
  }
  
  tt <- safe_t_test(a, b)
  pv <- if (is.null(tt)) NA_real_ else unname(tt$p.value)
  t_stat <- if (is.null(tt)) NA_real_ else unname(tt$statistic)
  sig_label <- p_to_label(pv)
  t_label <- if (is.na(t_stat)) "t = NA" else sprintf("t = %.2f", t_stat)
  
  y_max <- max(df$expression, na.rm = TRUE)
  if (!is.null(y_limit)) {
    y_bracket <- 1.82
    y_tick <- 1.68
    y_sig <- 2.22
    y_t <- 2.02
  } else {
    y_bracket <- if (y_max <= 0) 0.6 else y_max * 1.12
    y_tick <- if (y_max <= 0) 0.5 else y_max * 1.05
    y_sig <- y_bracket + max(0.18, y_max * 0.16)
    y_t <- y_bracket + max(0.06, y_max * 0.05)
  }
  sig_size <- if (identical(sig_label, "N/S")) 5 else 7
  
  p <- p +
    annotate("segment", x = x1, xend = x2, y = y_bracket, yend = y_bracket) +
    annotate("segment", x = x1, xend = x1, y = y_tick, yend = y_bracket) +
    annotate("segment", x = x2, xend = x2, y = y_tick, yend = y_bracket) +
    annotate(
      "text",
      x = mean(c(x1, x2)),
      y = y_sig,
      label = sig_label,
      vjust = 0.5,
      size = sig_size
    ) +
    annotate(
      "text",
      x = mean(c(x1, x2)),
      y = y_t,
      label = t_label,
      vjust = 0.5,
      size = 4.5
    )
  
  if (is.null(y_limit)) {
    p <- p + coord_cartesian(ylim = c(NA, y_sig * 1.12), clip = "off")
  } else {
    p <- p + coord_cartesian(ylim = c(-0.1, 2.45), clip = "off")
  }
  p
}

format_p <- function(p_value) {
  if (is.na(p_value)) {
    return("NA")
  }
  if (p_value < 0.001) {
    return(format(p_value, scientific = TRUE, digits = 3))
  }
  sprintf("%.4f", p_value)
}

build_anova_results <- function(df) {
  df <- df[!is.na(df$expression) & !is.na(df$group), ]
  df$group <- droplevels(factor(as.character(df$group)))
  groups <- levels(df$group)
  
  overall_p <- tryCatch(
    summary(aov(expression ~ group, data = df))[[1]]$`Pr(>F)`[1],
    error = function(e) NA_real_
  )
  
  combos <- combn(groups, 2, simplify = FALSE)
  pair_rows <- lapply(combos, function(pair) {
    pv <- safe_t_p(
      df$expression[df$group == pair[1]],
      df$expression[df$group == pair[2]]
    )
    data.frame(
      Comparison = paste(pair[1], "vs", pair[2]),
      `P-value` = format_p(pv),
      Significance = sig_status(pv),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })
  
  list(
    table = do.call(rbind, pair_rows),
    overall_p = overall_p,
    overall_label = sig_status(overall_p)
  )
}

make_anova_ui <- function(results) {
  rows <- lapply(seq_len(nrow(results$table)), function(i) {
    tags$tr(
      tags$td(results$table$Comparison[i]),
      tags$td(results$table[["P-value"]][i]),
      tags$td(results$table$Significance[i])
    )
  })
  
  tagList(
    tags$table(
      class = "table table-striped table-bordered",
      tags$thead(
        tags$tr(
          tags$th("Comparison"),
          tags$th("P-value"),
          tags$th("Significance")
        )
      ),
      tags$tbody(rows)
    ),
    tags$p(
      style = "margin-top: 12px; font-size: 1.15rem;",
      tags$strong("Overall significance: "),
      results$overall_label,
      " (p = ",
      format_p(results$overall_p),
      ")"
    )
  )
}

show_anova_modal <- function(title, results) {
  showModal(modalDialog(
    title = title,
    make_anova_ui(results),
    easyClose = TRUE,
    size = "m",
    footer = modalButton("Close")
  ))
}

make_vln <- function(obj, gene, fill_color, species_label, compare_groups = NULL) {
  
  validate(
    need(!is.na(gene), paste0("Undetectable expression, or no ", species_label, " homolog found."))
  )
  
  n_compare <- length(compare_groups)
  
  df <- fetch_gene_df(obj, gene)
  
  if (n_compare >= 2) {
    df <- df[as.character(df$group) %in% compare_groups, ]
    df$group <- droplevels(factor(as.character(df$group), levels = compare_groups))
  }
  
  p <- ggplot(df, aes(x = group, y = expression)) +
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
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 16),
      plot.margin = margin(12, 8, 8, 8)
    ) +
    labs(
      title = gene,
      x = NULL,
      y = "Lognorm Expression"
    )
  
  if (n_compare == 2) {
    p <- add_ttest_star(p, df, compare_groups[[1]], compare_groups[[2]])
  }
  
  p
}

wrap_across_label <- function(x) {
  gsub("^(Human: |Mouse: )", "\\1\n", x)
}

make_vln_from_df <- function(df, title) {
  df <- df[!is.na(df$expression) & !is.na(df$group), ]
  df$group <- droplevels(factor(df$group))
  
  validate(need(nrow(df) > 0, "No cells found for the selected categories."))
  
  n_groups <- nlevels(df$group)
  x_size <- if (n_groups >= 6) 11 else 12
  
  p <- ggplot(df, aes(x = group, y = expression, fill = species)) +
    geom_violin(
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
    scale_fill_manual(
      values = c(Human = human_color, Mouse = mouse_color),
      guide = "none"
    ) +
    scale_x_discrete(labels = wrap_across_label) +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 18),
      axis.text.x = element_text(size = x_size, lineheight = 0.95),
      axis.text.y = element_text(size = 16),
      plot.margin = margin(12, 16, 24, 8)
    ) +
    labs(
      title = title,
      x = NULL,
      y = "Lognorm Expression"
    )
  
  groups <- levels(df$group)
  p <- p +
    scale_y_continuous(
      breaks = seq(0, 2, by = 0.5)
    )
  
  if (length(groups) == 2) {
    p <- add_ttest_star(p, df, groups[[1]], groups[[2]], y_limit = 2)
  } else {
    p <- p + coord_cartesian(ylim = c(-0.1, 2), clip = "off")
  }
  
  p
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
  
  selected_source <- reactiveVal("source_Grinberg2026")
  
  observeEvent(input$source_Grinberg2026, {
    selected_source("source_Grinberg2026")
  })
  
  observeEvent(input$source_alzheimer, {
    selected_source("source_alzheimer")
  })
  
  observeEvent(input$source_schizophrenia, {
    selected_source("source_schizophrenia")
  })
  
  # Placeholder: Alzheimer and Schizophrenia reuse the Grinberg objects
  # until disease-specific datasets are added.
  current_human <- reactive({
    list(
      source_Grinberg2026 = human,
      source_alzheimer = human,
      source_schizophrenia = human
    )[[selected_source()]]
  })
  
  current_mouse <- reactive({
    list(
      source_Grinberg2026 = mouse,
      source_alzheimer = mouse,
      source_schizophrenia = mouse
    )[[selected_source()]]
  })
  
  cap_selection <- function(input_id, selected) {
    if (length(selected) > 8) {
      updateCheckboxGroupInput(session, input_id, selected = head(selected, 8))
      showNotification("Select at most 8 categories.", type = "warning")
      head(selected, 8)
    } else {
      selected
    }
  }
  
  observeEvent(input$human_clusters, {
    cap_selection("human_clusters", input$human_clusters)
  }, ignoreNULL = TRUE)
  
  observeEvent(input$mouse_clusters, {
    cap_selection("mouse_clusters", input$mouse_clusters)
  }, ignoreNULL = TRUE)
  
  across_active <- reactive({
    h <- input$human_clusters
    m <- input$mouse_clusters
    n_total <- length(h) + length(m)
    length(h) >= 1 && length(m) >= 1 && n_total >= 2 && n_total <= 8
  })
  
  across_dismissed <- reactiveVal(FALSE)
  
  dismiss_across <- function() {
    across_dismissed(TRUE)
    removeModal()
  }
  
  observeEvent(input$close_across_popup, {
    dismiss_across()
  })
  
  observeEvent(input$close_across_panel, {
    dismiss_across()
  })
  
  observeEvent(
    list(input$human_clusters, input$mouse_clusters, input$gene),
    {
      across_dismissed(FALSE)
      if (isTRUE(across_active()) && !is.null(input$gene) && nzchar(input$gene)) {
        n_cats <- length(input$human_clusters) + length(input$mouse_clusters)
        modal_w <- min(1800, 520 + n_cats * 170)
        plot_h <- if (n_cats >= 6) "400px" else "360px"
        showModal(modalDialog(
          title = "Across-category comparison",
          tags$style(HTML(sprintf(
            ".modal-dialog { max-width: min(96vw, %dpx) !important; width: min(96vw, %dpx) !important; }",
            modal_w,
            modal_w
          ))),
          plotOutput("across_vln_modal", height = plot_h),
          uiOutput("across_stats_modal"),
          easyClose = FALSE,
          size = "l",
          footer = actionButton(
            "close_across_popup",
            "Close",
            class = "btn btn-primary"
          )
        ))
      } else {
        removeModal()
      }
    },
    ignoreInit = TRUE
  )
  
  across_stats_ui <- function(n, anova_results) {
    if (n == 2) {
      return(star_note_ui())
    }
    if (n >= 3 && n <= 8 && !is.null(anova_results)) {
      return(div(class = "anova-panel", make_anova_ui(anova_results)))
    }
    NULL
  }
  
  across_combined_df <- reactive({
    req(isTRUE(across_active()))
    req(!is.null(input$gene), nzchar(input$gene))
    gp <- gene_pair()
    h_groups <- input$human_clusters
    m_groups <- input$mouse_clusters
    
    parts <- list()
    level_order <- character(0)
    
    if (length(h_groups) > 0 && !is.na(gp$human_gene)) {
      df_h <- fetch_gene_df(current_human(), gp$human_gene)
      df_h <- df_h[as.character(df_h$group) %in% h_groups, , drop = FALSE]
      df_h$group <- paste0("Human: ", as.character(df_h$group))
      df_h$species <- "Human"
      parts <- c(parts, list(df_h))
      level_order <- c(level_order, paste0("Human: ", h_groups))
    }
    
    if (length(m_groups) > 0 && !is.na(gp$mouse_gene)) {
      df_m <- fetch_gene_df(current_mouse(), gp$mouse_gene)
      df_m <- df_m[as.character(df_m$group) %in% m_groups, , drop = FALSE]
      df_m$group <- paste0("Mouse: ", as.character(df_m$group))
      df_m$species <- "Mouse"
      parts <- c(parts, list(df_m))
      level_order <- c(level_order, paste0("Mouse: ", m_groups))
    }
    
    validate(need(length(parts) > 0, "No matching categories found."))
    df <- do.call(rbind, parts)
    df$group <- factor(df$group, levels = unique(level_order))
    df
  })
  
  across_plot_title <- reactive({
    gp <- gene_pair()
    genes <- unique(c(gp$human_gene, gp$mouse_gene))
    genes <- genes[!is.na(genes)]
    paste(genes, collapse = " / ")
  })
  
  output$across_main_panel <- renderUI({
    if (
      !isTRUE(across_active()) ||
      isTRUE(across_dismissed()) ||
      is.null(input$gene) ||
      !nzchar(input$gene)
    ) {
      return(NULL)
    }
    n_cats_main <- length(input$human_clusters) + length(input$mouse_clusters)
    div(
      class = "across-main-panel",
      style = sprintf("min-width: %dpx; overflow-x: auto;", min(1600, 420 + n_cats_main * 150)),
      div(
        class = "across-main-header",
        h4("Across-category comparison"),
        actionButton(
          "close_across_panel",
          "Close",
          class = "btn btn-sm btn-outline-secondary"
        )
      ),
      plotOutput("across_vln", height = if (n_cats_main >= 6) "400px" else "360px"),
      uiOutput("across_stats")
    )
  })
  
  output$across_vln <- renderPlot({
    make_vln_from_df(across_combined_df(), across_plot_title())
  })
  
  output$across_vln_modal <- renderPlot({
    make_vln_from_df(across_combined_df(), across_plot_title())
  })
  
  output$across_stats <- renderUI({
    df <- across_combined_df()
    n <- nlevels(df$group)
    anova_res <- if (n >= 3 && n <= 8) build_anova_results(df) else NULL
    across_stats_ui(n, anova_res)
  })
  
  output$across_stats_modal <- renderUI({
    df <- across_combined_df()
    n <- nlevels(df$group)
    anova_res <- if (n >= 3 && n <= 8) build_anova_results(df) else NULL
    across_stats_ui(n, anova_res)
  })
  
  human_anova <- reactive({
    groups <- input$human_clusters
    if (length(groups) < 3 || length(groups) > 8) {
      return(NULL)
    }
    if (is.null(input$gene) || !nzchar(input$gene)) {
      return(NULL)
    }
    gp <- gene_pair()
    if (is.na(gp$human_gene)) {
      return(NULL)
    }
    df <- fetch_gene_df(current_human(), gp$human_gene)
    df <- df[as.character(df$group) %in% groups, ]
    build_anova_results(df)
  })
  
  mouse_anova <- reactive({
    groups <- input$mouse_clusters
    if (length(groups) < 3 || length(groups) > 8) {
      return(NULL)
    }
    if (is.null(input$gene) || !nzchar(input$gene)) {
      return(NULL)
    }
    gp <- gene_pair()
    if (is.na(gp$mouse_gene)) {
      return(NULL)
    }
    df <- fetch_gene_df(current_mouse(), gp$mouse_gene)
    df <- df[as.character(df$group) %in% groups, ]
    build_anova_results(df)
  })
  
  output$human_vln <- renderPlot({
    req(!is.null(input$gene),
        nzchar(input$gene),
        cancelOutput = TRUE
    )
    
    gp <- gene_pair()
    
    make_vln(
      obj = current_human(),
      gene = gp$human_gene,
      fill_color = human_color,
      species_label = "human",
      compare_groups = input$human_clusters
    )
  })
  
  output$mouse_vln <- renderPlot({
    req(!is.null(input$gene),
        nzchar(input$gene),
        cancelOutput = TRUE
    )
    
    gp <- gene_pair()
    
    make_vln(
      obj = current_mouse(),
      gene = gp$mouse_gene,
      fill_color = mouse_color,
      species_label = "mouse",
      compare_groups = input$mouse_clusters
    )
  })
  
  star_note_ui <- function() {
    p(
      class = "star-note",
      tags$strong("*"), " p < 0.05; ",
      tags$strong("**"), " p < 0.01; ",
      tags$strong("***"), " p < 0.001"
    )
  }
  
  output$human_star_note <- renderUI({
    if (length(input$human_clusters) != 2) {
      return(NULL)
    }
    star_note_ui()
  })
  
  output$mouse_star_note <- renderUI({
    if (length(input$mouse_clusters) != 2) {
      return(NULL)
    }
    star_note_ui()
  })
  
  output$human_anova_panel <- renderUI({
    results <- human_anova()
    if (is.null(results)) {
      return(NULL)
    }
    div(
      class = "anova-panel",
      make_anova_ui(results)
    )
  })
  
  output$mouse_anova_panel <- renderUI({
    results <- mouse_anova()
    if (is.null(results)) {
      return(NULL)
    }
    div(
      class = "anova-panel",
      make_anova_ui(results)
    )
  })
  
  # observeEvent(input$source_mouse_microglia, {
  #   selected_source("source_mouse_microglia")
  # })
  # 
  # observeEvent(input$source_mouse_astrocyte, {
  #   selected_source("source_mouse_astrocyte")
  # })
  
}