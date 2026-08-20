library(shiny)
library(bslib)

ui <- page_navbar(
  title = div(
    tags$img(
      src = "logo.png",
      style = "
      height: 120px;
      width: auto;
    "
    )
  ),
  
  window_title = "OPCEQ",
  
  theme = bs_theme(version = 5),
  
  header = tags$head(
    
    tags$style(HTML("
    
      .navbar {
        background-color: white !important;
        border-bottom: 1px solid #dee2e6;
        padding-top: 0px;
        padding-bottom: 0px;
      }
  
      .navbar .container-fluid {
        flex-direction: column;
        align-items: flex-start !important;
        gap: 0px;
      }
  
      .navbar-brand {
        margin-bottom: 0px;
        padding: 0;
        color: black !important;
      }
  
      .navbar-nav {
        width: 100%;
        margin-top: -20px;
  
      }
  
      .nav-link {
        color: black !important;
        font-size: 1.2rem;
      }
  
      .nav-link.active {
        border-bottom: none !important;
        box-shadow: none !important;
        font-weight: 700;
        color: #000 !important;
      }
      
      
      .selectize-control {
        width: 100% !important;
      }
  
      .selectize-input {
        width: 100% !important;
        min-height: 60px !important;
        font-size: 20px !important;
        padding: 12px 15px !important;
      }
      
      .selectize-input .item {
        font-size: 20px !important;
      }
  
      
      .source-list {
          margin-top: 15px;
        }
  
      .source-button {
        width: 100%;
        text-align: left;
        border: 2px solid transparent;
        border-radius: 0;
        background-color: #e9e9e9;
        color: #555;
        font-size: 1.1rem;
        font-weight: 600;
        padding: 15px;
        margin-bottom: 5px;
        white-space: normal;
      }
  
      .source-button:hover {
        background-color: white !important;
        color: #555 !important;
        border: 2px solid #777 !important;
      }
  
      .active-source {
        background-color: white !important;
        border: 2px solid #777 !important;
      }
  
      .source-button i {
        color: #2f75b5;
        margin-right: 10px;
      }
  
    "))
  ),
  
  nav_panel(
    "Gene Expression",
    
    fluidRow(
      
      column(
        width = 3,
        
        div(
          style = "
          position: sticky;
          top: 20px;
          height: 90vh;
          overflow-y: auto;
          padding-right: 15px;
        ",
          
          h4("Gene Search"),

          selectizeInput(
            "gene",
            "",
            choices = c(""),
            selected = "",
            options = list(
              placeholder = "GENE SYMBOL",
              create = FALSE,
              dropdownParent = "body",
              openOnFocus = FALSE,
              highlight = FALSE,
              maxOptions = 100,
              
              # Intercept the opening mechanism directly
              onDropdownOpen = I("
                function($dropdown) {
                  if (this.$control_input.val().length < 1) {
                    $dropdown.hide();
                  }
                }
              "),
              
              # Force drop down to hide/show dynamically as they backspace or type
              onType = I("
                function(str) {
                  if (str.length < 1) {
                    this.$dropdown.hide();
                  } else {
                    this.$dropdown.show();
                  }
                }
              ")
              
            )
          ),
          
          div(
            class = "source-list",
            
            actionButton(
              "source_Grinberg2026",
              label = tagList(icon("circle-chevron-right"), "Grinberg et al (2026)"),
              class = "source-button active-source"
            ),
            # 
            # actionButton(
            #   "source_mouse_microglia",
            #   label = tagList(icon("circle-chevron-right"), "Mouse Microglia Bennett et al. (2016)"),
            #   class = "source-button"
            # ),
            # 
            # actionButton(
            #   "source_mouse_astrocyte",
            #   label = tagList(icon("circle-chevron-right"), "Mouse Astrocyte Aging Clarke et al. (2018)"),
            #   class = "source-button"
            # ),
            
            tags$script(HTML("
                $(document).on('click', '.source-button', function() {
                  $('.source-button').removeClass('active-source');
                  $(this).addClass('active-source');
                });
              "))
          )
          
        )
 
      ),
      
      column(
        width = 9,
        
        card(
          full_screen = TRUE,
          card_header(
            tags$span(
              "Homo sapiens",
              style = "
                font-size: 1.5rem;
                font-weight: 500;
              "
            )
          ),
          plotOutput(
            "human_vln",
            height = "200px"
          )
        ),
        
        br(),
        
        card(
          full_screen = TRUE,
          card_header(
            tags$span(
              "Mus musculus",
              style = "
                font-size: 1.5rem;
                font-weight: 500;
              "
            )
          ),
          plotOutput(
            "mouse_vln",
            height = "200px"
          )
        )
      )
    )
  ),
  
  nav_panel(
    "Info",
    # h3("Dataset information"),
    p("The correspondence between human OPC clusters and the labels defined by Grinberg et al. (2026) is shown: "),
    tags$table(
      class = "table table-striped table-bordered",
      
      tags$thead(
        tags$tr(
          tags$th("This website"),
          tags$th("Grinberg et al. (2026) Figure 1")
        )
      ),
      
      tags$tbody(
        tags$tr(
          tags$td("Developmental OPC"),
          tags$td("eCluster 1")
        ),
        tags$tr(
          tags$td("Proliferative OPC (G1-to-S phase)"),
          tags$td("eCluster 3")
        ),
        tags$tr(
          tags$td("Proliferative OPC"),
          tags$td("eCluster 2")
        ),
        tags$tr(
          tags$td("Synaptic OPC"),
          tags$td("eCluster 4")
        )
      )
    ),
    
    br(),
    
    p("The correspondence between the mouse OPC clusters in OPCeq and the clusters defined in the data sources is as follows:"),
    
    tags$table(
      class = "table table-striped table-bordered",
      
      tags$thead(
        tags$tr(
          tags$th("This website"),
          tags$th("Data sources")
        )
      ),
      
      tags$tbody(
        tags$tr(
          tags$td("Active OPC"),
          tags$td("Dennis et al. (2024)")
        ),
        tags$tr(
          tags$td("Proliferative OPC"),
          tags$td("Dennis et al. (2024); Heo et al. (2025)")
        ),
        tags$tr(
          tags$td("Homeostatic OPC"),
          tags$td("Dennis et al. (2024)")
        )
      )
    )
  ),
  
  nav_panel(
    "Download",
    # h3("Downloads"),
    p("The raw FASTQ files are available via the NCBI Gene Expression Omnibus (GEO) under sample ID: GSM6598840.")
  ),
  
  nav_panel(
    "Citation",
    # h3("Citation"),
    p("Please cite Grinberg et al. within any publication that makes use of analyses inspired by OPCeq."),
    p("Grinberg AR, Liu L, Patel I, Guntur B, Glass IA; Birth Defects Research Laboratory; Shinoda K, Nobuta H. Transcriptional Heterogeneity Reveals a Synaptic Gene Program in Developing and Adult Human Oligodendrocyte Precursor Cells. bioRxiv [Preprint]. 2026 Jan 1:2026.01.01.697017. doi: 10.64898/2026.01.01.697017. PMID: 41509247; PMCID: PMC12776387."),
    p("Dennis DJ, Wang BS, Karamboulas K, Kaplan DR, Miller FD. Single-cell approaches define two groups of mammalian oligodendrocyte precursor cells and their evolution over developmental time. Stem Cell Reports. 2024 May 14;19(5):654-672. doi: 10.1016/j.stemcr.2024.03.002. Epub 2024 Apr 4. PMID: 38579710; PMCID: PMC11103788."),
    p("Heo D, Kim AA, Neumann B, Doze VN, Xu YKT, Mironova YA, Slosberg J, Goff LA, Franklin RJM, Bergles DE. Transcriptional profiles of mouse oligodendrocyte precursor cells across the lifespan. Nat Aging. 2025 Apr;5(4):675-690. doi: 10.1038/s43587-025-00840-2. Epub 2025 Mar 31. PMID: 40164771.")
  ),
  
  nav_panel(
    "Contact",
    # h3("Contact"),
    p("For general enquiries or comments about the website, please contact:"),
    tags$p(
      tags$strong("Hiroko Nobuta"),
      br(),
      tags$a(
        href = "mailto:nobuta@cabm.rutgers.edu",
        "nobuta@cabm.rutgers.edu"
      ),
      br(),
      tags$a(
        href = "https://cabm.rutgers.edu/research/nobuta-lab",
        target = "_blank",
        "https://cabm.rutgers.edu/research/nobuta-lab"
      )
    ),
    p("For dataset specific questions, please contact:"),
    tags$p(
      tags$strong("Kosaku Shinoda"),
      br(),
      tags$a(
        href = "mailto:kosaku.shinoda@einsteinmed.edu",
        "kosaku.shinoda@einsteinmed.edu"
      )
    )
    
  )
)