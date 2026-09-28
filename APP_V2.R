#README:

#Hello, this app was created to facilitate investigation into fungal antimicrovial resistance
#worldwide. The project started by analyzing the FungAMR database which still has its own database
#we then decided to build another database to compare with the FungAMR database before deciding making
#this app public and allowing others to add data would help build an interactive tool
#for those of us who are researching fungal antimicrobial resistance. This is my first
#coding project and was done iteratively over several session heavily relying on stack overflow and googling how to
# do various things. Leading to several reused blocks. ChatGPT was used to test for bughunting and stress testing. 
#This involved the AI finding errors in the code and suggesting fixes, which did not work majority of the time, 
# which were then reviewed by myself and use as a reference when fixing the bugs 
#Claude also provided an example of how to allow users to upload data, which was used as a reference
#in a similar manner, meaning that while google, stackoverflow, and AI agents were used as references and as places to ask questions the coding was done by me.
#This also means that there are some silly errors and areas of code that can probably be removed safely, however the app is in a working
#state so this is not a priority currently. Eventually the code will be cleaned and streamlined in a later release.
#Thank you for using the app, I hope you learn something from the data and contribute some of your own.

#-------------------------------------------------------------------------------
#Download Dependencies (For local run)
#-------------------------------------------------------------------------------
#install.packages("shiny")
#install.packages("bslib")
#install.packages("dplyr")
#install.packages("ggplot2")
#install.packages("rnaturalearth")
#install.packages("sf")
#install.packages("stringr")
#install.packages("leaflet")
#install.packages("purrr")
#install.packages("DT")
#install.packages("ggalluvial")
#install.packages("readxl")
#install.packages("viridis")
#install.packages("tidyr")
#install.packages("rsconnect")
#install.packages("RSQLite")
#install.packages("shinyjs")
#-------------------------------------------------------------------------------
#Load Dependencies (For Local Run)
#-------------------------------------------------------------------------------
library(shiny)
library(bslib)
library(dplyr)
library(ggplot2)
library(rnaturalearth)
library(sf)
library(stringr)
library(leaflet)
library(purrr)
library(DT)
library(ggalluvial)
library(readxl)
library(viridis)
library(tidyr)
library(rsconnect)
library(RSQLite)
library(shinyjs)

#-------------------------------------------------------------------------------
#Internal database
#-------------------------------------------------------------------------------

db <- dbConnect(RSQLite::SQLite(),"fungal_amr.db")

if(!dbExistsTable(db,"submissions")) {
  submissions_schema <- data.frame(
    'Family/species' = character(),
    Gene = character(),
    Mutation = character(),
    Source_Type = character(),
    Location = character(),
    country = character(),
    Accession = character(),
    sequence = character(),
    submission_status = character(),
    check.names = FALSE
  )
  dbCreateTable(db, "submissions", submissions_schema)
} else {
  existing_cols <- dbListFields(db, "submissions")
  review_cols <- c("submission_status")
  for (col in setdiff(review_cols, existing_cols)) {
    dbExecute(db, paste0('ALTER TABLE submissions ADD COLUMN "', col, '" TEXT'))
  }
  
}

#-------------------------------------------------------------------------------
#Validation
#-------------------------------------------------------------------------------

validate_and_normalize_upload <- function(df) {
  messages <- c()
  required_cols <- c(
    "Family/species","Gene","Mutation","Source_Type","Location","Accession","sequence"
  )
  
  
  if(!all(required_cols %in% colnames(df)) ||
     any(df[required_cols] == "" | is.na(df[required_cols])) ||
     nrow(df) == 0) { return("You are Missing Data")}
  
  df_normalized <- df
  
  df_normalized$country <- df_normalized$Location
  df_normalized$submission_status <- "pending_review"
  
  
  return(df_normalized) }

#-------------------------------------------------------------------------------
# Load Data
#-------------------------------------------------------------------------------

Fungal_AMR <- read.csv("Data/Composite_FungAMR_ONLYGEODATA.csv", header = TRUE)

Taxonomy <- read.csv("Data/fungal_taxonomy.csv",header = TRUE)

Fungal_AMR_Taxonomy <- Fungal_AMR %>%
  mutate(species = trimws(tolower(species))) %>%
  left_join(
    Taxonomy %>%
      mutate(species = trimws(tolower(species))),
    by = "species"
  )

Comp_Fungal_AMR <- read_xlsx("Data/Manually_Curated_Dataset.xlsx")

Public_Data <- read_xlsx("Data/Public_Data.xlsx")

#-------------------------------------------------------------------------------
# Data Normalization - FungAMR
#-------------------------------------------------------------------------------
Fungal_AMR_Taxonomy <- Fungal_AMR_Taxonomy %>%
  mutate(
    mutation = str_trim(mutation),
    mutation = case_when(
      mutation %in% c("F129L,F129L", "F129L") ~ "F129L",
      TRUE ~ mutation
    ),
    gene.or.protein.name = str_trim(
      tolower(gene.or.protein.name)
    ),
    gene.or.protein.name = case_when(
      gene.or.protein.name %in% c(
        "sqle",
        "squalene epoxidase"
      ) ~ "sqle",
      gene.or.protein.name %in% c(
        "cytochrome b",
        "cytb"
      ) ~ "cytb",
      gene.or.protein.name %in% c(
        "cyp51a"
      ) ~ "cyp51a",
      gene.or.protein.name %in% c(
        "cyp51"
      ) ~ "cyp51",
      gene.or.protein.name %in% c(
        "erg11"
      ) ~ "erg11",
      gene.or.protein.name %in% c(
        "tub2",
        "beta-tubulin 2"
      ) ~ "beta-tubulin 2",
      TRUE ~ gene.or.protein.name
    ))
#-------------------------------------------------------------------------------
# Data Prep - FungAMR
#-------------------------------------------------------------------------------
Fungal_AMR_Taxonomy <- Fungal_AMR_Taxonomy %>%
  mutate(
    country = if_else(
      is.na(geographic_region) | geographic_region == "",
      NA_character_,
      str_split(geographic_region, ":") %>% map_chr(1)
    ),
    country = case_when(
      country == "USA"                   ~ "United States of America",
      country == "Brazil Brazil"         ~ "Brazil",
      country == "China,Nanjing,Jiangsu" ~ "China",
      country == "UK"                    ~ "United Kingdom",
      country == "China "                ~ "China",
      TRUE                               ~ country
    )
  )
FungAMR_Clinical <- Fungal_AMR_Taxonomy %>%
  filter(strain.origin.if.available == "Clinical")

FungAMR_Agricultural <- Fungal_AMR_Taxonomy %>%
  filter(strain.origin.if.available == "Environment")

FungAMR_Artificiall <- Fungal_AMR_Taxonomy %>%
  filter(strain.origin.if.available == "Lab")

FungAMR_Misc <- Fungal_AMR_Taxonomy %>%
  filter(strain.origin.if.available %in% c("Unknown", "Evolved"))



#-------------------------------------------------------------------------------
# Data Normalization - FEL
#-------------------------------------------------------------------------------
Comp_Fungal_AMR <- Comp_Fungal_AMR %>%
  mutate(
    Gene = str_trim(tolower(Gene)),
    Gene = case_when(
      Gene %in% c(
        "cyp51a"
      ) ~ "cyp51a",
      Gene %in% c(
        "cytochrome b",
        "cytb"
      ) ~ "cytb",
      Gene %in% c(
        "erg11"
      ) ~ "erg11",
      Gene %in% c(
        "tub2",
        "beta-tubulin 2"
      ) ~ "beta-tubulin 2",
      Gene %in% c(
        "cyp51"
      ) ~ "cyp51",
      TRUE ~ Gene
    ))

Comp_Fungal_AMR <- Comp_Fungal_AMR %>%
  mutate( `Family/species` = case_when (`Family/species` %in%
                                          c("Alternaria sp") ~ "Pleosporaceae",
                                        `Family/species` %in%
                                          c("Aspergillus", "Aspergillus fumigatus") ~ "Aspergillaceae",
                                        `Family/species` %in%
                                          c("Nectriacieae") ~ "Nectriaceae",
                                        `Family/species` %in%
                                          c("Colletotrichum truncatum") ~ "Glomerellaceae",
                                        `Family/species` %in%
                                          c("Lasiodiplodia theobromae") ~ "Botryosphaeriaceae",
                                        `Family/species` %in%
                                          c("Pseudocercospora fijiensis") ~ "Mycosphaerellaceae",
                                        TRUE ~ `Family/species`))
#-------------------------------------------------------------------------------
#Data Normalization Public Data
#-------------------------------------------------------------------------------

Public_Data <- Public_Data %>%
  mutate(
    country =  if_else( is.na(Location), 
                        NA_character_, 
                        Location ))



#-------------------------------------------------------------------------------
# DATA Prep - FEL
#-------------------------------------------------------------------------------
Comp_Fungal_AMR <- Comp_Fungal_AMR %>%
  mutate(
    country = if_else(
      is.na(Location) | Location == "",
      NA_character_,
      str_split(Location, ":") %>% map_chr(1)
    ),
    country = case_when(
      country == "USA"                   ~ "United States of America",
      country == "Brazil Brazil"         ~ "Brazil",
      country == "China,Nanjing,Jiangsu" ~ "China",
      country == "UK"                    ~ "United Kingdom",
      country == "China "                ~ "China",
      country == "japan"                 ~ "Japan",
      country == "japan "                ~ "Japan",
      TRUE                               ~ country
    ))

FEL_Clinical <- Comp_Fungal_AMR %>%
  filter(Source_Type == "Clinical")

FEL_Agricultural <- Comp_Fungal_AMR %>%
  filter(Source_Type == "Agricultural")

FEL_Animal <- Comp_Fungal_AMR %>%
  filter(Source_Type == "Animal")

FEL_Soil <- Comp_Fungal_AMR %>%
  filter(Source_Type == "Soil")

FEL_Artificial <- Comp_Fungal_AMR %>%
  filter(Source_Type == "Artificial")

FEL_Combined <- Comp_Fungal_AMR

#-------------------------------------------------------------------------------
#Public Source Data
#-------------------------------------------------------------------------------
Public_Data_Clinical <- Public_Data %>%
  filter(Source_Type == "Clinical")

Public_Data_Agricultural <- Public_Data %>%
  filter(Source_Type == "Agricultural")

Public_Data_Misc <- Public_Data %>%
  filter(Source_Type == "Misc")
#-------------------------------------------------------------------------------
# UI
#-------------------------------------------------------------------------------
ui <- page_sidebar(
  useShinyjs(),
  title = tags$span(
    tags$img(
      src = "FELLogo.png",
      height = "40px",
      style = "margin-right: 15px; vertical-align: middle;"
    ),
    tags$span(
      "FungAMR & FEL Database Comparison",
      style = "vertical-align: middle; font-size: 24px;"
    )),
  theme = bs_theme(
    bg = "#FFFFFF",
    fg = "#000000",
    primary = "#0d6efd",
    base_font = font_google("Inter")
  ),
  sidebar = sidebar(
    helpText("Data by Category"),
    radioButtons(
      inputId = "page_selection",
      label = NULL,
      choices = c(
        "Home"                    = "home",
        "How to Add Data"         =  "instructions",
        "Methods"                 = "methods",
        "Clinical Data"           = "clinical",
        "Agricultural Data"       = "agricultural",
        "Resistance Genes"        = "genes",
        "Mutations"               = "mutations",
        "Geographic Distribution" = "geographic",
        "Data by Family"          = "family",
        "Compare Databases"       = "compare"
      ),
      selected = "home"
    )),
  uiOutput("page_content"))
#-------------------------------------------------------------------------------
# Server
#-------------------------------------------------------------------------------
options(shiny.maxRequestSize = 50*1024^2)

server <- function(input, output, session) {
  
 
  is_admin <- reactiveVal(FALSE)
  refresh_trigger <- reactiveVal(0)
  

  active_page <- reactiveVal("home")
  
  observeEvent(input$page_selection, {
    active_page(input$page_selection)
  })
  
  observeEvent(input$goto_review, {
    active_page("review")
  })
  
  observeEvent(input$back_to_home, {
    active_page("home")
  })
  
  approved_user_submissions <- reactive({
    refresh_trigger()
    tryCatch({
      df <- dbGetQuery(db, "SELECT * FROM submissions WHERE submission_status = 'approved' ")
      if (nrow(df) == 0) {
        return(data.frame())
      }
      df %>%
        select(`Family/species`, Gene, Mutation, Source_Type, Location, country)
    }, error = function(e) {
      data.frame()
    })
  })
  
  public_data_combined <- reactive({
    bind_rows(Public_Data, approved_user_submissions())
  })
  
  public_clinical_r <- reactive({
    public_data_combined() %>% filter(Source_Type == "Clinical")
  })
  
  public_agricultural_r <- reactive({
    public_data_combined() %>% filter(Source_Type == "Agricultural")
  })
  
  public_misc_r <- reactive({
    public_data_combined() %>% filter(Source_Type == "Misc")
  })
  
  combined_geo_data <- reactive({
    bind_rows(
      Fungal_AMR_Taxonomy %>% select(country, 'Family/species' = Family),
      FEL_Combined %>% select(country, 'Family/species'),
      public_data_combined() %>% select(country, 'Family/species')
    )
  })
  
  pending_submissions <- reactive({
    refresh_trigger()
    tryCatch({
      dbGetQuery(db, "SELECT rowid AS submission_id, * FROM submissions WHERE submission_status = 'pending_review'")
    }, error = function(e) {
      data.frame()
    })
  })
  
  reviewed_submissions <- reactive({
    refresh_trigger()
    tryCatch({
      dbGetQuery(db, "SELECT rowid AS submission_id, * FROM submissions WHERE submission_status != 'pending_review'")
    }, error = function(e) {
      data.frame()
    })
  })
  
  make_count_table <- function(df, col, label, page_length = 10) {
    tryCatch({
      result <- df %>%
        filter(!is.na(.data[[col]])) %>%
        count(!!sym(label) := .data[[col]], sort = TRUE)
      
      datatable(
        result,
        options = list(pageLength = page_length)
      )
    },
    error = function(e) {
      datatable(
        data.frame(
          Error = paste("Could not load data:", e$message)
        ))})}
  
  make_geographic_map <- function(df, title) {
    country_counts <- df %>%
      filter(!is.na(country)) %>%
      count(country, name = "samples")
    
    world <- ne_countries(
      scale = "medium",
      returnclass = "sf")
    
    world_data <- left_join(
      world,
      country_counts,
      by = c("name" = "country"))
    
    ggplot(world_data) +
      geom_sf(
        aes(fill = samples),
        color = "grey90") +
      scale_fill_viridis(
        option = "C",
        na.value = "white"
      ) +
      theme_minimal() +
      labs(
        title = title,
        fill = "Number of Samples")}
  
  
  make_interactive_geographic_map <- function(
    df,
    family_column,
    title)
  {country_data <- df %>%
    filter(!is.na(country), country != "") %>%
    group_by(country) %>%
    summarise(
      families = paste(
        sort(unique(.data[[family_column]][
          !is.na(.data[[family_column]]) &
            .data[[family_column]] != ""
        ])),
        collapse = "<br>"
      ),
      family_count = n_distinct(
        .data[[family_column]],
        na.rm = TRUE
      ),
      .groups = "drop"
    )
  
  world <- ne_countries(
    scale = "medium",
    returnclass = "sf"
  )
  
  world_data <- world %>%
    left_join(
      country_data,
      by = c("name" = "country")
    )
  
  world_data <- world_data %>%
    mutate(
      families_display = if_else(families %in% c(NA, ""), "No Families Recorded",
                                 families)
      ,
      
      family_count_display = case_when(
        is.na(family_count) ~ "0",
        TRUE ~ as.character(family_count)
      ),
      
      hover_text = paste0(
        "<strong>",
        name,
        "</strong>",
        "<br><br>",
        "<strong>Families:</strong> ",
        family_count_display,
        "<br><br>",
        "<strong>Families found:</strong><br>",
        families_display
      )
    )
  
  leaflet(
    world_data,
    options = leafletOptions(
      minZoom = 1,
      maxZoom = 8
    )
  ) %>%
    addProviderTiles(
      providers$CartoDB.Positron
    ) %>%
    addPolygons(
      fillColor = ~ifelse(
        is.na(family_count),
        "#BDBDBD",
        "#08519C"
      ),
      fillOpacity = 0.75,
      color = "white",
      weight = 0.8,
      opacity = 1,
      label = ~lapply(
        hover_text,
        HTML
      ),
      highlightOptions = highlightOptions(
        weight = 2,
        color = "#000000",
        fillOpacity = 0.9,
        bringToFront = TRUE
      )
    ) %>%
    addControl(
      html = paste0(
        "<div style='font-size:18px;font-weight:bold;'>",
        title,
        "</div>"
      ),
      position = "topright"
    ) %>%
    addLegend(
      position = "bottomright",
      colors = c("#08519c", "#BDBDBD"),
      labels = c("Data present", "No Data"),
      opacity = 2
    )}
  
  uploaded_data_validated <- reactiveVal(NULL)
  
  output$download_template <- downloadHandler(
    filename = "FungAMR_Upload_Template.csv",
    content = function(file) {
      template <- data.frame(
        `Family/species` = c("Aspergillaceae", "Nectriaceae", "Pleosporaceae"),
        Gene = c("cyp51a", "cyp51", "cyp51a"),
        Mutation = c("G448S", "Y136F", "F129L"),
        Source_Type = c("Clinical", "Agricultural", "Animal"),
        Location = c("United States of America", "China", "Brazil"),
        Accession = c("", "", ""),
        sequence = c("", "", ""),
        check.names = FALSE
      )
      write.csv(template, file, row.names = FALSE)
    }
  )
  
  observeEvent(input$submit_upload, {
    req(input$user_upload)
    
    tryCatch({
      raw_data <- read.csv(input$user_upload$datapath, stringsAsFactors = FALSE, check.names = FALSE)
      validation_result <- validate_and_normalize_upload(raw_data)
      
      if (is.character(validation_result)) {
        output$validation_results <- renderUI({
          div(class = "alert alert-danger", h4("validation Failed"), p(validation_result))
        })
        output$upload_preview <- renderDT(NULL)
        return()
      }
      
      uploaded_data_validated(validation_result)
      
      output$validation_results <- renderUI({
        div(
          class = "alert alert-info",
          h4("Data Validated Successfully"),
          p(strong(paste("Rows to be added:", nrow(validation_result)))),
          br(),
      
          actionButton("confirm_upload", "Add to Database", class = "btn-success btn-lg"),
          actionButton("cancel_upload", "Cancel", class = "btn-secondary btn-lg")
        )
      })
      
      output$upload_preview <- renderDT({
        preview_data <- validation_result %>%
          select(`Family/species`, Gene, Mutation, Source_Type, Location, country)
        datatable(preview_data, caption = "Preview of normalized data (first 10 rows)", options = list(pageLength = 10))
      })
      
    }, error = function(e) {
      output$validation_results <- renderUI({
        div(class = "alert alert-danger", h4("Error Reading File"), p(paste("Error:", e$message)))
      })
      output$upload_preview <- renderDT(NULL)
    })
  })
  
  observeEvent(input$confirm_upload, {
    req(uploaded_data_validated())
    
    
      validated_data <- uploaded_data_validated()
      dbAppendTable(db, "submissions", validated_data)
      refresh_trigger(refresh_trigger() + 1)
      
      output$validation_results <- renderUI({
        div(
          class = "alert alert-success",
          h4("Success!"),
          p(strong(nrow(validated_data)), "records have been submitted to the database."),
          p("Your submission is pending review. Thank you for contributing!")
        )
      })
      
      output$upload_preview <- renderDT(NULL)
      reset("user_upload")
      uploaded_data_validated(NULL)
  })
  
  observeEvent(input$cancel_upload, {
    uploaded_data_validated(NULL)
    output$validation_results <- renderUI(NULL)
    output$upload_preview <- renderDT(NULL)
    reset("user_upload")
  })
  

  
  output$review_panel <- renderUI({
    if (!is_admin()) {
      div(
        h2("Review Submissions"),
        p("This section is restricted to trusted reviewers. Enter the reviewer password to continue."),
        passwordInput("admin_password", "Reviewer Password:"),
        actionButton("admin_login", "Log In", class = "btn-primary"),
        actionButton("back_to_home", "Back to Home", class = "btn-link"),
        br(), br(),
        uiOutput("admin_login_message")
      )
    } else {
      div(
        h2("Review Submissions"),
        p("Approve or reject data submitted by users. Approved submissions are added to the Public Data used throughout the app."),
        actionButton("admin_logout", "Log Out", class = "btn-secondary"),
        actionButton("back_to_home", "Back to Home", class = "btn-link"),
        br(), br(),
        h4("Pending Submissions"),
        DTOutput("pending_submissions_table"),
        br(),
        actionButton("approve_selected", "Approve Selected", class = "btn-success"),
        actionButton("reject_selected", "Reject Selected", class = "btn-danger"),
        br(), br(),
        h4("Review History"),
        DTOutput("reviewed_submissions_table")
      )
    }
  })
  
  observeEvent(input$admin_login, {
    reviewer_password <- Sys.getenv("REVIEWER_PASSWORD")
    if (!is.null(input$admin_password) && nzchar(input$admin_password) && input$admin_password == reviewer_password) {
      is_admin(TRUE)
      output$admin_login_message <- renderUI(NULL)
    } else {
      output$admin_login_message <- renderUI({
        div(class = "alert alert-danger", "Incorrect password.")
      })
    }
  })
  
  observeEvent(input$admin_logout, {
    is_admin(FALSE)
  })
  
  output$pending_submissions_table <- renderDT({
    datatable(
      pending_submissions(),
      rownames = FALSE,
      selection = "multiple",
      options = list(pageLength = 10, scrollX = TRUE)
    )
  })
  
  output$reviewed_submissions_table <- renderDT({
    datatable(
      reviewed_submissions(),
      rownames = FALSE,
      options = list(pageLength = 10, scrollX = TRUE)
    )
  })
  
  observeEvent(input$approve_selected, {
    req(is_admin())
    selected_rows <- input$pending_submissions_table_rows_selected
    req(length(selected_rows) > 0)
    
    ids <- pending_submissions()$submission_id[selected_rows]
    for (id in ids) {
      dbExecute(db, "UPDATE submissions SET submission_status = 'approved' WHERE rowid = ?", params = list(id))
    }
    refresh_trigger(refresh_trigger() + 1)
    showNotification(paste(length(ids), "submission(s) approved."), type = "message")
  })
  
  observeEvent(input$reject_selected, {
    req(is_admin())
    selected_rows <- input$pending_submissions_table_rows_selected
    req(length(selected_rows) > 0)
    
    ids <- pending_submissions()$submission_id[selected_rows]
    for (id in ids) {
      dbExecute(db, "UPDATE submissions SET submission_status = 'rejected' WHERE rowid = ?", params = list(id))
    }
    refresh_trigger(refresh_trigger() + 1)
    showNotification(paste(length(ids), "submission(s) rejected."), type = "warning")
  })
  
  output$page_content <- renderUI({
    switch(
      active_page(),
      
      
      
      
      "home" = div(
        h2("Welcome to the Database Comparison Tool"),
        p("This Shiny app allows you to explore and compare FungAMR and FEL databases with the abillity to provide additional data to a public dataset for comparisons."),
        p("Select a category from the sidebar to begin."),
        br(),
        h3("Data Loading Status"),
        
        if (nrow(Fungal_AMR_Taxonomy) > 0) {
          p(
            style = "color:green",
            " FungAMR data loaded successfully"
          )
        } else {
          p(
            style = "color:red",
            " Unsuccsesful in Loading FungAMR data "
          )
        },
        
        if (nrow(Comp_Fungal_AMR) > 0) {
          p(
            style = "color:green",
            "FEL data loaded successfully"
          )
        } else {
          p(
            style = "color:red",
            "Unsuccesful in loading FungAMR data"
          )
        },
        
        if (nrow(public_data_combined()) > 0) {
          p(
            style = "color:green",
            " Public data loaded successfully"
          )
        } else {
          p(
            style = "color:red",
            " Unsuccsesful in Loading Public data "
          )
        },
        
        
        br(),
        h3("Database Overview"),
        
        p(
          "FungAMR Database: ",
          nrow(Fungal_AMR_Taxonomy),
          " records"
        ),
        
        p(
          "  - Clinical: ",
          nrow(FungAMR_Clinical)
        ),
        
        p(
          "  - Agricultural: ",
          nrow(FungAMR_Agricultural)
        ),
        
        p(
          "  - Artifical: ",
          nrow(FungAMR_Artificiall)
        ),
        br(),
        
        p(
          "  - Misc: ",
          nrow(FungAMR_Misc)
        ),
        br(),
        
        p(
          "FEL Database: ",
          nrow(Comp_Fungal_AMR),
          " records"
        ),
        
        p(
          "  - Clinical: ",
          nrow(FEL_Clinical)
        ),
        
        p(
          "  - Agricultural: ",
          nrow(FEL_Agricultural)
        ),
        
        p(
          "  - Animal: ",
          nrow(FEL_Animal)
        ),
        
        p(
          "  - Soil: ",
          nrow(FEL_Soil)
        ),
        
        p(
          "  - Artificial: ",
          nrow(FEL_Artificial)
        ),
        br(),
        
        p(
          "Public Data:",
          nrow(public_data_combined())
        ),
        
        p(
          " - Artificial",
          nrow(public_agricultural_r())
        ),
        
        p(
          "  - Clinical",
          nrow(public_clinical_r())
        ),
        
        p(
          "  - Misc",
          nrow(public_misc_r())
        ),
        
        tags$div(
          style = "margin-top: 60px;",
          tags$small(
            actionLink(
              "goto_review",
              "Reviewer access",
              style = "color:#aaaaaa; font-size:11px; text-decoration:underline;"
            )
          )
        ),
      ),
      
      "instructions" = div(
        h2("Add Your Data to Our Database"),
        p("We welcome contributions! Please follow these steps:"),
        br(),
        h4("Step 1: Download the Template"),
        p("Start with our template to ensure your data is formatted correctly."),
        downloadButton(
          "download_template",
          "Download Template CSV",
          class = "btn-primary btn-lg"
        ),
        br(), br(),
        h4("Step 2: Fill in Your Data"),
        p("Required columns:"),
        tags$ul(
          tags$li(strong("Family/species"), " - Fungal family or species name"),
          tags$li(strong("Gene"), " - Resistance gene (e.g., cyp51a, cytb, erg11)"),
          tags$li(strong("Mutation"), " - Specific mutation (e.g., G448S, F129L)"),
          tags$li(strong("Source_Type"), " - Clinical, Agricultural, Animal, Soil, Artificial, or Misc"),
          tags$li(strong("Location"), " - Country or geographic location"),
          tags$li(strong("Accession"), " - Provide an ncbi accession number"),
          tags$li(strong("Sequence"), " - Please provide the associated sequence containing the mutations")
        ),
        br(),
        h4("Step 3: Upload Your File"),
        fileInput(
          inputId = "user_upload",
          label = "Select CSV file",
          accept = c(".csv"),
          buttonLabel = "Browse..."
        ),
        actionButton(
          inputId = "submit_upload",
          label = "Validate and Preview Data",
          class = "btn-primary btn-lg"
        ),
        br(), br(),
        uiOutput("validation_results"),
        DTOutput("upload_preview")
      ),
      
      
      "methods" = div(
        h1("  Methodology"),
        br(),
        h3(
          style = "color:blue",
          "Database"),
        p("This project started after reading the FungAMR paper by Bedard et al. published in Nature on the 11th of April 2025. The goal
        was just to map the data they provided in their database to discover if there were location specific and global patterns,
        of fungal antimicrobial resistance. To do this we utilized the accessions provided in their database and extracted
        geographic data that was attached to the accessions. This resulted in a signficant loss of data in comparison to the original database.
        so we decided to build another database to act as a comparison tool and to help identify more patterns. Once we start to explore this we decided to make the
        tool we were using public to encourage others to supply data to an already fully interactive tool."),

        br(),
        h3(
          style = "color:blue", "Shiny App"
        ),
        p("To build this database the original codes I was utilizing to examine global fungicide resistance data in R
        were modified for use in a Shiny Application. They were then made interactive to make it easier for users unfamiliar with the
        data to explore. All coding was done by hand but referenced stack overflow, several websites, ChatGPT and Claude for help with structuring
        and when utilizing new R plugins. A github stores all publically available version of the app and allows for people to make suggestions or keave comments."),
      ),

      "clinical" = div(
        h2("Clinical Data"),
        
        tabsetPanel(
          
          tabPanel(
            "FungAMR Overview",
            p(
              "FungAMR Clinical samples: ",
              nrow(FungAMR_Clinical)
            ),
            DTOutput("fungamr_clinical_table")
          ),
          
          tabPanel(
            "FEL Overview",
            p(
              "FEL Clinical samples: ",
              nrow(FEL_Clinical)
            ),
            DTOutput("fel_clinical_table")
          ),
          
          tabPanel(
            "Publi Data Overview",
            p(
              "Public Data Clinical samples: ",
              nrow(FEL_Clinical)
            ),
            DTOutput("public_clinical_table")
          ),
          
          tabPanel(
            "Geographic Map - FungAMR",
            plotOutput(
              "fungamr_clinical_map",
              height = "600px"
            )
          ),
          
          tabPanel(
            "Geographic Map - FEL",
            plotOutput(
              "fel_clinical_map",
              height = "600px"
            )
          ),
          tabPanel(
            "Geographic Map - Public",
            plotOutput(
              "Public_clinical_map",
              height = "600px"
            )
          ),
          
          tabPanel(
            "Genes Comparison",
            p("FungAMR Clinical Genes:"),
            DTOutput("fungamr_clinical_genes"),
            br(),
            p("FEL Clinical Genes:"),
            DTOutput("fel_clinical_genes"),
            br(),
            p("Public Data Genes"),
            DTOutput("Public_Data_Clinical_genes")
          ),
          
          tabPanel(
            "Mutations Comparison",
            p("FungAMR Clinical Mutations:"),
            DTOutput("fungamr_clinical_mutations"),
            br(),
            p("FEL Clinical Mutations:"),
            DTOutput("fel_clinical_mutations"),
            br(),
            p("Public Data Genes"),
            DTOutput("Public_Data_Clinical_mutations")
          )
        )
      ),
      
      
      "agricultural" = div(
        h2("Agricultural Data"),
        
        tabsetPanel(
          
          tabPanel(
            "FungAMR Overview",
            p(
              "FungAMR Agricultural samples: ",
              nrow(FungAMR_Agricultural)
            ),
            DTOutput("fungamr_agricultural_table")
          ),
          
          tabPanel(
            "Public Data Overview",
            p(
              "Public Data Agricultural Samples: ",
              nrow(public_agricultural_r())
            ),
            DTOutput("Public_Data_Agricultural_table")
          ),
          
          tabPanel(
            "FEL Overview",
            p(
              "FEL Agricultural samples: ",
              nrow(FEL_Agricultural)
            ),
            DTOutput("fel_agricultural_table")
          ),
          
          tabPanel(
            "Geographic Map - FungAMR",
            plotOutput(
              "fungamr_agricultural_map",
              height = "600px"
            )
          ),
          
          tabPanel(
            "Geographic Map - FEL",
            plotOutput(
              "fel_agricultural_map",
              height = "600px"
            )
          ),
          
          tabPanel(
            "Geographic Map - Public Data",
            plotOutput(
              "PublicData_Agricultural_Map"
            ),
          ),
          
          tabPanel(
            "Genes Comparison",
            p("FungAMR Agricultural Genes:"),
            DTOutput("fungamr_agricultural_genes"),
            br(),
            p("FEL Agricultural Genes:"),
            DTOutput("fel_agricultural_genes"),
            br(),
            p("Public Data Agricultural Genes"),
            DTOutput("Public_Data_Agricultural_Genes")
          ),
          
          tabPanel(
            "Mutations Comparison",
            p("FungAMR Agricultural Mutations:"),
            DTOutput("fungamr_agricultural_mutations"),
            br(),
            p("FEL Agricultural Mutations:"),
            DTOutput("fel_agricultural_mutations"),
            br(),
            p("Public Data Agricultural Mutations:"),
            DTOutput("Public_Data_Agricultural_mutations")
          )
        )
      ),
      
      
      
      "genes" = div(
        h2("Resistance Genes"),
        
        tabsetPanel(
          
          tabPanel(
            "FungAMR Data",
            
            h4("FungAMR Resistance Genes"),
            p(
              "Resistance genes identified in the FungAMR database."
            ),
            
            DTOutput("fungamr_resistance_genes")
          ),
          
          tabPanel(
            "FEL Data",
            
            h4("FEL Resistance Genes"),
            p(
              "Resistance genes identified in the FEL database."
            ),
            
            DTOutput("fel_resistance_genes")
          ),
          
          tabPanel(
            "Public Data",
            
            h4("Public Data Resistance Genes"),
            p(
              "Resistance genes identified in the publically build database."
            ),
            
            DTOutput("Public_Data_resistance_genes")
          ),
          
          tabPanel(
            "Comparison",
            
            h4("FungAMR vs FEL VS Public Data Resistance Genes"),
            p(
              "Comparison of resistance genes identified in each database."
            ),
            
            DTOutput("resistance_gene_comparison")
          )
        )
      ),
      
      
      "mutations" = div(
        h2("Mutations"),
        
        tabsetPanel(
          
          tabPanel(
            "FungAMR Data",
            
            h4("FungAMR Resistance Mutations"),
            p(
              "Resistance-associated mutations identified in the FungAMR database."
            ),
            
            DTOutput("fungamr_resistance_mutations")
          ),
          
          tabPanel(
            "FEL Data",
            
            h4("FEL Resistance Mutations"),
            p(
              "Resistance-associated mutations identified in the FEL database."
            ),
            
            DTOutput("fel_resistance_mutations")
          ),
          
          tabPanel(
            "Public Data",
            
            h4("Public Data Resistance Mutations"),
            p(
              "Resistance-associated mutations identified in the Public database."
            ),
            
            DTOutput("Public_Data_resistance_mutations")
          ),
          
          tabPanel(
            "Comparison",
            
            h4("FungAMR vs FEL Vs Public Data Resistance Mutations"),
            p(
              "Comparison of resistance-associated mutations identified in each database."
            ),
            
            DTOutput("resistance_mutation_comparison")
          )
        )
      ),
      
      
      
      "geographic" = div(
        h2("Geographic Distribution(Please Give Time To Load)"),
        
        radioButtons(
          "geo_type",
          "Select Origin Type:",
          choices = c(
            "Clinical" = "clinical",
            "Agricultural" = "agricultural"
          ),
          inline = TRUE
        ),
        
        br(),
        
        tabsetPanel(
          
          tabPanel(
            "FungAMR Map",
            leafletOutput(
              "fungamr_geo_map",
              height = "600px"
            )
          ),
          
          
          
          tabPanel(
            "FEL Map",
            leafletOutput(
              "fel_geo_map",
              height = "600px"
            )
          ),
          
          tabPanel(
            "public Data Map",
            leafletOutput(
              "Public_Data_geo_map",
              height = "600px"
            )
          ),

          tabPanel(
            "Combined Data map",
            leafletOutput(
              "combined_geo_map",
              height = "600px"
            )
          )
          
        )
      ),
      
      
      
      
      
      "family" = div(
        h2("Data by Family"),
        
        uiOutput("family_select"),
        
        br(),
        
        tabsetPanel(
          
          tabPanel(
            "FungAMR Alluvial",
            plotOutput(
              "fungamr_alluvial",
              height = "600px"
            )
          ),
          
          tabPanel(
            "FEL Alluvial",
            plotOutput(
              "fel_alluvial",
              height = "600px"
            )
          ),
          tabPanel(
            "Public Data Alluvial",
            plotOutput(
              "Public_Data_alluvial",
              height = "600px"
            )
          )
        )
      ),
      
      
      
      "compare" = div(
        h2("Compare Databases"),
        
        tabsetPanel(
          
          tabPanel(
            "Summary Statistics",
            
            h4("Database Size Comparison"),
            
            p(
              "FungAMR Total Records: ",
              nrow(Fungal_AMR_Taxonomy)
            ),
            
            p(
              "FEL Total Records: ",
              nrow(Comp_Fungal_AMR)
            ),
            
            p(
              "Public Data Records:",
              nrow(public_data_combined())
            ),
            
            br(),
            
            h4("Clinical Records"),
            
            p(
              "FungAMR Clinical: ",
              nrow(FungAMR_Clinical)
            ),
            
            p(
              "FEL Clinical: ",
              nrow(FEL_Clinical)
            ),
            
            p("Public Data Clinical:",
              nrow(public_clinical_r())
              
            ),
            
            br(),
            
            h4("Agricultural Records"),
            
            p(
              "FungAMR Agricultural: ",
              nrow(FungAMR_Agricultural)
            ),
            
            p(
              "FEL Agricultural: ",
              nrow(FEL_Agricultural)
            ),
            
            p(
              "Public Data Clinical:",
              nrow(public_agricultural_r())
            ),
            
          ),
          
          tabPanel(
            "Geographic Comparison",
            
            h4("FungAMR Geographic Distribution"),
            
            plotOutput(
              "fungamr_compare_map",
              height = "500px"
            ),
            
            br(),
            
            h4("FEL Geographic Distribution"),
            
            plotOutput(
              "fel_compare_map",
              height = "500px"
            )
          ),
          
          tabPanel(
            "Family Distribution",
            
            h4("FungAMR Top Families"),
            
            plotOutput(
              "fungamr_family_dist",
              height = "500px"
            ),
            
            br(),
            
            h4("FEL Top Families"),
            
            plotOutput(
              "fel_family_dist",
              height = "500px"
            )
          ),
          
          tabPanel(
            "Resistance Genes",
            
            h4("Gene Comparison - FungAMR vs FEL"),
            
            DTOutput("compare_gene_table")
          ),
          
          tabPanel(
            "Mutations",
            
            h4("Mutation Comparison - FungAMR vs FEL"),
            
            DTOutput("compare_mutation_table")
          )
        )
      ),
      
      "review" = uiOutput("review_panel")
    )
  })
  
  
  output$fungamr_clinical_table <- renderDT(
    make_count_table(
      FungAMR_Clinical,
      "Family",
      "Family",
      
    )
  )
  
  output$fel_clinical_table <- renderDT(
    make_count_table(
      FEL_Clinical,
      "Family/species",
      "Family",
      
    )
  )
  
  output$public_clinical_table <- renderDT(
    make_count_table(
      public_clinical_r(),
      "Family/species",
      "Family",
      
    )
  )
  
  output$fungamr_clinical_genes <- renderDT(
    make_count_table(
      FungAMR_Clinical,
      "gene.or.protein.name",
      "Gene",
      
    )
  )
  
  output$fel_clinical_genes <- renderDT(
    make_count_table(
      FEL_Clinical,
      "Gene",
      "Gene",
      
    )
  )
  
  output$fungamr_clinical_mutations <- renderDT(
    make_count_table(
      FungAMR_Clinical,
      "mutation",
      "Mutation"
    )
  )
  
  output$fel_clinical_mutations <- renderDT(
    make_count_table(
      FEL_Clinical,
      "Mutation",
      "Mutation"
    )
  )
  
  output$Public_Data_Clinical_genes <- renderDT(
    make_count_table(
      public_clinical_r(),
      "Gene",
      "Gene"
    )
  )
  
  
  
  output$Public_Data_Clinical_mutations <- renderDT(
    make_count_table(
      public_clinical_r(),
      "Mutation",
      "Mutation"
    )
  )
  
  output$fungamr_clinical_map <- renderPlot({
    make_geographic_map(
      FungAMR_Clinical,
      "FungAMR Clinical Geographic Distribution"
    )
  })
  
  output$fel_clinical_map <- renderPlot({
    make_geographic_map(
      FEL_Clinical,
      "FEL Clinical Geographic Distribution"
    )
  })
  
  output$Public_clinical_map <- renderPlot({
    make_geographic_map(
      public_clinical_r(),
      "Public Data Clinical Geographic Distribution"
    )
  })
  
  
  output$fungamr_agricultural_table <- renderDT(
    make_count_table(
      FungAMR_Agricultural,
      "Family",
      "Family"
    )
  )
  
  output$fel_agricultural_table <- renderDT(
    make_count_table(
      FEL_Agricultural,
      "Family/species",
      "Family"
    )
  )
  
  output$fungamr_agricultural_genes <- renderDT(
    make_count_table(
      FungAMR_Agricultural,
      "gene.or.protein.name",
      "Gene"
    )
  )
  
  output$Public_Data_Agricultural_table <- renderDT(
    make_count_table(
      public_agricultural_r(),
      "Family/species",
      "Family"
    )
  )
  
  output$fel_agricultural_genes <- renderDT(
    make_count_table(
      FEL_Agricultural,
      "Gene",
      "Gene"
    )
  )
  
  output$Public_Data_Agricultural_Genes <- renderDT(
    make_count_table(
      public_agricultural_r(),
      "Gene",
      "Gene"
    )
  )
  
  output$fungamr_agricultural_mutations <- renderDT(
    make_count_table(
      FungAMR_Agricultural,
      "mutation",
      "Mutation"
    )
  )
  
  output$fel_agricultural_mutations <- renderDT(
    make_count_table(
      FEL_Agricultural,
      "Mutation",
      "Mutation"
    )
  )
  
  output$Public_Data_Agricultural_mutations <- renderDT(
    make_count_table(
      public_agricultural_r(),
      "Mutation",
      "Mutation"
    )
  )
  
  output$fungamr_agricultural_map <- renderPlot({
    make_geographic_map(
      FungAMR_Agricultural,
      "FungAMR Agricultural Geographic Distribution"
    )
  })
  
  output$fel_agricultural_map <- renderPlot({
    make_geographic_map(
      FEL_Agricultural,
      "FEL Agricultural Geographic Distribution"
    )
  })
  
  output$PublicData_Agricultural_Map <- renderPlot({
    make_geographic_map(
      public_agricultural_r(),
      "Public Data Agricultural Geographic Distribution"
    )
  })
  
  
  output$fungamr_resistance_genes <- renderDT({
    
    gene_data <- Fungal_AMR_Taxonomy %>%
      filter(
        !is.na(gene.or.protein.name),
        gene.or.protein.name != ""
      ) %>%
      group_by(gene.or.protein.name) %>%
      summarise(
        Records = n(),
        Species = n_distinct(species, na.rm = TRUE),
        Families = n_distinct(Family, na.rm = TRUE),
        Countries = n_distinct(country, na.rm = TRUE),
        Mutations = n_distinct(
          mutation[
            !is.na(mutation) &
              mutation != ""
          ]
        ),
        .groups = "drop"
      ) %>%
      rename(
        Gene = gene.or.protein.name
      ) %>%
      arrange(desc(Records))
    
    datatable(
      gene_data,
      rownames = FALSE,
      options = list(
        pageLength = 15,
        scrollX = TRUE
      )
    )
  })
  
  
  output$fel_resistance_genes <- renderDT({
    
    gene_data <- FEL_Combined %>%
      filter(
        !is.na(Gene),
        Gene != ""
      ) %>%
      group_by(Gene) %>%
      summarise(
        Records = n(),
        Families_or_Species = n_distinct(
          `Family/species`,
          na.rm = TRUE
        ),
        Countries = n_distinct(
          country,
          na.rm = TRUE
        ),
        Mutations = n_distinct(
          Mutation[
            !is.na(Mutation) &
              Mutation != ""
          ]
        ),
        Sources = n_distinct(
          Source_Type,
          na.rm = TRUE
        ),
        .groups = "drop"
      ) %>%
      arrange(desc(Records))
    
    datatable(
      gene_data,
      rownames = FALSE,
      options = list(
        pageLength = 15,
        scrollX = TRUE
      )
    )
  })
  
  output$Public_Data_resistance_genes <- renderDT({
    
    gene_data <- public_data_combined() %>%
      filter(
        !is.na(Gene),
        Gene != ""
      ) %>%
      group_by(Gene) %>%
      summarise(
        Records = n(),
        Families_or_Species = n_distinct(
          `Family/species`,
          na.rm = TRUE
        ),
        Countries = n_distinct(
          country,
          na.rm = TRUE
        ),
        Mutations = n_distinct(
          Mutation[
            !is.na(Mutation) &
              Mutation != ""
          ]
        ),
        Sources = n_distinct(
          Source_Type,
          na.rm = TRUE
        ),
        .groups = "drop"
      ) %>%
      arrange(desc(Records))
    
    datatable(
      gene_data,
      rownames = FALSE,
      options = list(
        pageLength = 15,
        scrollX = TRUE
      )
    )
  })
  
  output$resistance_gene_comparison <- renderDT({
    
    fungamr_genes <- Fungal_AMR_Taxonomy %>%
      filter(
        !is.na(gene.or.protein.name),
        gene.or.protein.name != ""
      ) %>%
      count(
        gene.or.protein.name,
        name = "FungAMR_Count"
      ) %>%
      rename(
        Gene = gene.or.protein.name
      )
    
    fel_genes <- FEL_Combined %>%
      filter(
        !is.na(Gene),
        Gene != ""
      ) %>%
      count(
        Gene,
        name = "FEL_Count"
      )
    
    Public_Data_Genes <- public_data_combined() %>%
      filter(
        !is.na(Gene),
        Gene != ""
      ) %>%
      count(
        Gene,
        name = "Public_Data_Count"
      )
    
    comparison <- full_join(
      fungamr_genes,
      fel_genes,
      by = "Gene"
    ) %>%
      full_join(
        Public_Data_Genes,
        by = "Gene"
      )%>%
      mutate(
        FungAMR_Count = replace_na(
          FungAMR_Count,
          0
        ),
        FEL_Count = replace_na(
          FEL_Count,
          0
        ),
        Public_Data_Count = replace_na(
          Public_Data_Count,
          0
        ),
        
        Total = FungAMR_Count + FEL_Count+Public_Data_Count,
        Database = case_when(
          FungAMR_Count > 0 &
            Public_Data_Count > 0 &
            FEL_Count > 0 ~ "All",
          FungAMR_Count > 0 ~ "FungAMR Only",
          Public_Data_Count >0 ~ "Public Data",
          FEL_Count > 0 ~ "FEL Only",
          TRUE ~ "None"
        )
      ) %>%
      arrange(desc(Total))
    
    datatable(
      comparison,
      rownames = FALSE,
      options = list(
        pageLength = 15,
        scrollX = TRUE
      )
    )
  })
  
  
  
  output$fungamr_resistance_mutations <- renderDT({
    
    mutation_data <- Fungal_AMR_Taxonomy %>%
      filter(
        !is.na(mutation),
        mutation != ""
      ) %>%
      group_by(mutation) %>%
      summarise(
        Records = n(),
        Genes = n_distinct(
          gene.or.protein.name[
            !is.na(gene.or.protein.name) &
              gene.or.protein.name != ""
          ]
        ),
        Species = n_distinct(
          species,
          na.rm = TRUE
        ),
        Families = n_distinct(
          Family,
          na.rm = TRUE
        ),
        Countries = n_distinct(
          country,
          na.rm = TRUE
        ),
        .groups = "drop"
      ) %>%
      rename(
        Mutation = mutation
      ) %>%
      arrange(desc(Records))
    
    datatable(
      mutation_data,
      rownames = FALSE,
      options = list(
        pageLength = 15,
        scrollX = TRUE
      )
    )
  })
  
  
  output$fel_resistance_mutations <- renderDT({
    
    mutation_data <- FEL_Combined %>%
      filter(
        !is.na(Mutation),
        Mutation != ""
      ) %>%
      group_by(Mutation) %>%
      summarise(
        Records = n(),
        Genes = n_distinct(
          Gene[
            !is.na(Gene) &
              Gene != ""
          ]
        ),
        Families_or_Species = n_distinct(
          `Family/species`,
          na.rm = TRUE
        ),
        Countries = n_distinct(
          country,
          na.rm = TRUE
        ),
        Sources = n_distinct(
          Source_Type,
          na.rm = TRUE
        ),
        .groups = "drop"
      ) %>%
      arrange(desc(Records))
    
    datatable(
      mutation_data,
      rownames = FALSE,
      options = list(
        pageLength = 15,
        scrollX = TRUE
      )
    )
  })
  
  output$Public_Data_resistance_mutations <- renderDT({
    
    mutation_data <- public_data_combined() %>%
      filter(
        !is.na(Mutation),
        Mutation != ""
      ) %>%
      group_by(Mutation) %>%
      summarise(
        Records = n(),
        Genes = n_distinct(
          Gene[
            !is.na(Gene) &
              Gene != ""
          ]
        ),
        Families_or_Species = n_distinct(
          `Family/species`,
          na.rm = TRUE
        ),
        Countries = n_distinct(
          country,
          na.rm = TRUE
        ),
        Sources = n_distinct(
          Source_Type,
          na.rm = TRUE
        ),
        .groups = "drop"
      ) %>%
      arrange(desc(Records))
    
    datatable(
      mutation_data,
      rownames = FALSE,
      options = list(
        pageLength = 15,
        scrollX = TRUE
      )
    )
  })
  
  output$resistance_mutation_comparison <- renderDT({
    
    fungamr_mutations <- Fungal_AMR_Taxonomy %>%
      filter(
        !is.na(mutation),
        mutation != ""
      ) %>%
      count(
        mutation,
        name = "FungAMR_Count"
      ) %>%
      rename(
        Mutation = mutation
      )
    
    fel_mutations <- FEL_Combined %>%
      filter(
        !is.na(Mutation),
        Mutation != ""
      ) %>%
      count(
        Mutation,
        name = "FEL_Count"
      )
    
    Public_Mutations <- public_data_combined() %>%
      filter( !is.na(Mutation),
              Mutation != ""
      ) %>%
      count(
        Mutation,
        name = "Public_Data_Count"
      )
    
    comparison <- full_join(
      fungamr_mutations,
      fel_mutations,
      by = "Mutation"
    ) %>%
      full_join(
        Public_Mutations,
        by = "Mutation"
      ) %>%
      mutate(
        FungAMR_Count = replace_na(
          FungAMR_Count,
          0
        ),
        FEL_Count = replace_na(
          FEL_Count,
          0
          
        ),
        Public_Data_Count = replace_na(
          Public_Data_Count,
          0
        ),
        
        Total = FungAMR_Count + FEL_Count+ Public_Data_Count,
        Database = case_when(
          FungAMR_Count > 0 &
            Public_Data_Count > 0 &
            FEL_Count > 0 ~ "All",
          FungAMR_Count > 0 ~ "FungAMR Only",
          Public_Data_Count > 0 ~ "Public Data Only",
          FEL_Count > 0 ~ "FEL Only",
          TRUE ~ "None"
        )
      ) %>%
      arrange(desc(Total))
    
    datatable(
      comparison,
      rownames = FALSE,
      options = list(
        pageLength = 15,
        scrollX = TRUE
      )
    )
  })
  
  
  
  output$fungamr_geo_map <- renderLeaflet({
    
    if (input$geo_type == "clinical") {
      
      make_interactive_geographic_map(
        FungAMR_Clinical,
        "Family",
        "FungAMR Clinical Geographic Distribution"
      )
      
    } else {
      
      make_interactive_geographic_map(
        FungAMR_Agricultural,
        "Family",
        "FungAMR Agricultural Geographic Distribution"
      )
    }
  })
  
  
  output$fel_geo_map <- renderLeaflet({
    
    if (input$geo_type == "clinical") {
      
      make_interactive_geographic_map(
        FEL_Clinical,
        "Family/species",
        "FEL Clinical Geographic Distribution"
      )
      
    } else {
      
      make_interactive_geographic_map(
        FEL_Agricultural,
        "Family/species",
        "FEL Agricultural Geographic Distribution"
      )
    }
  })
  
  output$Public_Data_geo_map <- renderLeaflet({
    
    if (input$geo_type == "clinical") {
      
      make_interactive_geographic_map(
        public_clinical_r(),
        "Family/species",
        "Public Data Clinical Geographic Distribution"
      )
      
    } else {
      
      make_interactive_geographic_map(
        public_agricultural_r(),
        "Family/species",
        "Public Data Agricultural Geographic Distribution"
      )
    }
  })

  output$combined_geo_map <- renderLeaflet({
    make_interactive_geographic_map(
      combined_geo_data(),
      "Family/species",
      "Combined Geographic Distribution (All Databases)"
    )
  })


  output$family_select <- renderUI({
    
    fungamr_families <- Fungal_AMR_Taxonomy %>%
      filter(!is.na(Family)) %>%
      pull(Family) %>%
      unique() %>%
      sort()
    
    fel_families <- FEL_Combined %>%
      filter(!is.na(`Family/species`)) %>%
      pull(`Family/species`) %>%
      unique() %>%
      sort()
    
    public_families <- public_data_combined() %>%
      filter(!is.na(`Family/species`)) %>%
      pull(`Family/species`) %>%
      unique() %>%
      sort()
    
    all_families <- c(
      fungamr_families,
      fel_families,
      public_families
    ) %>%
      unique() %>%
      sort()
    
    selectInput(
      "selected_family",
      "Select a Family:",
      choices = all_families
    )
  })
  
  
  output$fungamr_alluvial <- renderPlot({
    
    req(input$selected_family)
    
    data_plot <- Fungal_AMR_Taxonomy %>%
      filter(
        Family == input$selected_family,
        !is.na(country)
      ) %>%
      group_by(
        Family,
        country,
        gene.or.protein.name,
        mutation
      ) %>%
      summarise(
        count = n(),
        .groups = "drop"
      )
    
    if (nrow(data_plot) == 0) {
      
      return(
        ggplot() +
          annotate(
            "text",
            x = 0.5,
            y = 0.5,
            label = "No data available for this family",
            size = 5
          ) +
          theme_void()
      )
    }
    
    ggplot(
      data_plot,
      aes(
        axis1 = Family,
        axis2 = country,
        axis3 = gene.or.protein.name,
        axis4 = mutation,
        y = count
      )
    ) +
      geom_alluvium(
        aes(fill = gene.or.protein.name),
        alpha = 0.8
      ) +
      geom_stratum() +
      geom_text(
        stat = "stratum",
        aes(label = after_stat(stratum)),
        size = 3
      ) +
      scale_x_discrete(
        limits = c(
          "Family",
          "Country",
          "Gene",
          "Mutation"
        ),
        expand = c(.1, .1)
      ) +
      scale_y_continuous(
        expand = c(0, 0)
      ) +
      theme_bw() +
      labs(
        title = paste(
          "FungAMR -",
          input$selected_family
        ),
        y = "Count",
        fill = "Gene"
      )
  })
  
  
  output$fel_alluvial <- renderPlot({
    
    req(input$selected_family)
    
    data_plot <- FEL_Combined %>%
      filter(
        `Family/species` == input$selected_family,
        !is.na(country)
      ) %>%
      group_by(
        `Family/species`,
        country,
        Gene,
        Mutation
      ) %>%
      summarise(
        count = n(),
        .groups = "drop"
      )
    
    if (nrow(data_plot) == 0) {
      
      return(
        ggplot() +
          annotate(
            "text",
            x = 0.5,
            y = 0.5,
            label = "No data available for this family",
            size = 5
          ) +
          theme_void()
      )
    }
    
    ggplot(
      data_plot,
      aes(
        axis1 = `Family/species`,
        axis2 = country,
        axis3 = Gene,
        axis4 = Mutation,
        y = count
      )
    ) +
      geom_alluvium(
        aes(fill = Gene),
        alpha = 0.8
      ) +
      geom_stratum() +
      geom_text(
        stat = "stratum",
        aes(label = after_stat(stratum)),
        size = 3
      ) +
      scale_x_discrete(
        limits = c(
          "Family",
          "Country",
          "Gene",
          "Mutation"
        ),
        expand = c(.1, .1)
      ) +
      scale_y_continuous(
        expand = c(0, 0)
      ) +
      theme_bw() +
      labs(
        title = paste(
          "FEL -",
          input$selected_family
        ),
        y = "Count",
        fill = "Gene"
      )
  })
  
  
  output$Public_Data_alluvial <- renderPlot({
    
    req(input$selected_family)
    
    data_plot <- public_data_combined() %>%
      filter(
        `Family/species` == input$selected_family,
        !is.na(country)
      ) %>%
      group_by(
        `Family/species`,
        country,
        Gene,
        Mutation
      ) %>%
      summarise(
        count = n(),
        .groups = "drop"
      )
    
    if (nrow(data_plot) == 0) {
      
      return(
        ggplot() +
          annotate(
            "text",
            x = 0.5,
            y = 0.5,
            label = "No data available for this family",
            size = 5
          ) +
          theme_void()
      )
    }
    
    ggplot(
      data_plot,
      aes(
        axis1 = `Family/species`,
        axis2 = country,
        axis3 = Gene,
        axis4 = Mutation,
        y = count
      )
    ) +
      geom_alluvium(
        aes(fill = Gene),
        alpha = 0.8
      ) +
      geom_stratum() +
      geom_text(
        stat = "stratum",
        aes(label = after_stat(stratum)),
        size = 3
      ) +
      scale_x_discrete(
        limits = c(
          "Family",
          "Country",
          "Gene",
          "Mutation"
        ),
        expand = c(.1, .1)
      ) +
      scale_y_continuous(
        expand = c(0, 0)
      ) +
      theme_bw() +
      labs(
        title = paste(
          "Public -",
          input$selected_family
        ),
        y = "Count",
        fill = "Gene"
      )
  })
  
  output$fungamr_compare_map <- renderPlot({
    
    make_geographic_map(
      Fungal_AMR_Taxonomy,
      "FungAMR Geographic Distribution"
    )
  })
  
  
  output$fel_compare_map <- renderPlot({
    
    make_geographic_map(
      FEL_Combined,
      "FEL Geographic Distribution"
    )
  })
  
  
  output$fungamr_family_dist <- renderPlot({
    
    family_counts <- Fungal_AMR_Taxonomy %>%
      filter(!is.na(Family)) %>%
      count(
        Family,
        sort = TRUE
      ) %>%
      arrange(desc(n)) %>%
      slice(1:15)
    
    ggplot(
      family_counts,
      aes(
        x = reorder(Family, n),
        y = n
      )
    ) +
      geom_col(
        fill = "#08519c"
      ) +
      coord_flip() +
      labs(
        title = "FungAMR - Top 15 Families",
        x = "Family",
        y = "Count"
      ) +
      theme_minimal()
  })
  
  
  output$fel_family_dist <- renderPlot({
    
    family_counts <- FEL_Combined %>%
      filter(!is.na(`Family/species`)) %>%
      count(
        `Family/species`,
        sort = TRUE
      ) %>%
      arrange(desc(n)) %>%
      slice(1:15)
    
    ggplot(
      family_counts,
      aes(
        x = reorder(`Family/species`, n),
        y = n
      )
    ) +
      geom_col(
        fill = "#d8b365"
      ) +
      coord_flip() +
      labs(
        title = "FEL - Top 15 Families",
        x = "Family",
        y = "Count"
      ) +
      theme_minimal()
  })
  
  
  output$compare_gene_table <- renderDT({
    
    fungamr_genes <- Fungal_AMR_Taxonomy %>%
      filter(!is.na(gene.or.protein.name)) %>%
      count(
        gene.or.protein.name,
        sort = TRUE
      ) %>%
      rename(
        Gene = gene.or.protein.name,
        FungAMR_Count = n
      )
    
    fel_genes <- FEL_Combined %>%
      filter(!is.na(Gene)) %>%
      count(
        Gene,
        sort = TRUE
      ) %>%
      rename(
        FEL_Count = n
      )
    
    Public_Genes <- public_data_combined() %>%
      filter(!is.na(Gene)) %>%
      count(
        Gene,
        sort = TRUE
      ) %>%
      rename(
        Public_Count = n
      )
    
    comparison <- full_join(
      fungamr_genes,
      fel_genes,
      by = "Gene"
    ) %>%
      full_join(
        Public_Genes,
        by = "Gene"
      )%>%
      mutate(
        FungAMR_Count = replace_na(
          FungAMR_Count,
          0
        ),
        FEL_Count = replace_na(
          FEL_Count,
          0
        ),
        Public_Count = replace_na (
          Public_Count,
          0
        )) %>%
      arrange(desc(FungAMR_Count))
    
    datatable(
      comparison,
      options = list(
        pageLength = 15
      )
    )
  })
  
  
  output$compare_mutation_table <- renderDT({
    
    fungamr_mutations <- Fungal_AMR_Taxonomy %>%
      filter(!is.na(mutation)) %>%
      count(
        mutation,
        sort = TRUE
      ) %>%
      rename(
        Mutation_Name = mutation,
        FungAMR_Count = n
      )
    
    fel_mutations <- FEL_Combined %>%
      filter(!is.na(Mutation)) %>%
      count(
        Mutation,
        sort = TRUE
      ) %>%
      rename(
        Mutation_Name = Mutation,
        FEL_Count = n
      )
    
    Public_mutations <- public_data_combined() %>%
      filter(!is.na(Mutation)) %>%
      count(
        Mutation,
        sort = TRUE
      ) %>%
      rename(
        Mutation_Name = Mutation,
        Public_Count = n
      )
    
    comparison<- full_join(
      fungamr_mutations,
      fel_mutations,
      by = "Mutation_Name"
    ) %>%
      full_join(
        Public_mutations,
        by = "Mutation_Name"
      ) %>%
      mutate(
        FungAMR_Count = replace_na(
          FungAMR_Count,
          0
        ),
        FEL_Count = replace_na(
          FEL_Count,
          0
        ),
        Public_Count = replace_na(
          Public_Count,
          0
        ),
        
        Difference = FungAMR_Count - FEL_Count
      ) %>%
      arrange(desc(FungAMR_Count))
    
    datatable(
      comparison,
      options = list(
        pageLength = 15
      )
    )
  })
}

onStop(function()  { dbDisconnect(db)})

#-------------------------------------------------------------------------------
# Run
#-------------------------------------------------------------------------------

shinyApp(
  ui,
  server
)
