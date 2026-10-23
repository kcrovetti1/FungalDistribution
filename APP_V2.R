#-------------------------------------------------------------------------------
#Necessary Libraries - set up
#-------------------------------------------------------------------------------

### This section lists all of the libraries that are needed to run this shiny app
### Install them using install.packages("package name") if needed

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
library(tidyr)
library(RSQLite)
library(shinyjs)
library(bsicons)
library(ggiraph)
library(mapgl)

#-------------------------------------------------------------------------------
#building a combined dataset - set up
#-------------------------------------------------------------------------------
Combined_Data <- read_xlsx ("Data/Combined_Database.xlsx")



Taxonomic_data <- read.csv("Data/fungal_taxonomy.csv",header = TRUE)

Combined_Data <- Combined_Data %>%
  mutate(species = trimws(tolower(`Family/species`))) %>%
  left_join(
    Taxonomic_data %>%
      mutate(species = trimws(tolower(species))),
    by = "species") %>%
  mutate(Family = coalesce(Family, `Family/species` ))



Combined_Taxonomy<- Combined_Data %>%
  mutate(
    Mutation = str_trim(Mutation),
    Mutation = case_when(
      Mutation %in% c("F129L,F129L", "F129L") ~ "F129L",
      TRUE ~ Mutation
    ),
    Gene = str_trim(
      tolower(Gene)
    ),
    Gene = case_when(
      Gene %in% c(
        "sqle",
        "squalene epoxidase"
      ) ~ "sqle",
      Gene %in% c(
        "cytochrome b",
        "cytb"
      ) ~ "cytb",
      Gene %in% c(
        "cyp51a"
      ) ~ "cyp51a",
      Gene %in% c(
        "erg11"
      ) ~ "erg11",
      Gene %in% c(
        "tub2",
        "beta-tubulin 2"
      ) ~ "beta-tubulin 2",
      Gene %in% c(
        "b tubulin",
        "beta-tubulin"
      ) ~ "beta-tubulin",
      Gene %in% c(
        "ffcyp51b",
        "cyp51b"
      ) ~ "cyp51b",
      Gene %in% c(
        "ltcyp51",
        "cyp51"
      ) ~ "cyp51",
      TRUE ~ Gene
    ))




Combined_Taxonomy<- Combined_Taxonomy%>%
  mutate( Family = case_when (Family %in%
                                c("Alternaria sp") ~ "Pleosporaceae",
                              Family %in%
                                c("Aspergillus", "Aspergillus fumigatus") ~ "Aspergillaceae",
                              Family %in%
                                c("Nectriacieae") ~ "Nectriaceae",
                              Family%in%
                                c("Colletotrichum truncatum") ~ "Glomerellaceae",
                              Family %in%
                                c("Lasiodiplodia theobromae") ~ "Botryosphaeriaceae",
                              Family %in%
                                c("Pseudocercospora fijiensis") ~ "Mycosphaerellaceae",
                              TRUE ~ Family))

Combined_Taxonomy<- Combined_Taxonomy %>%
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
      country == "brazil"                ~ "Brazil",
      country == "Isreal"                ~ "Israel",
      country == "Suadi Arabia"          ~"Saudi Arabia" ,
      country =="Czech Republic" ~ "Czechia",
      TRUE                               ~ country
    ))


#-------------------------------------------------------------------------------
#Removing duplicated and set up images - set up
#-------------------------------------------------------------------------------

### This step removed duplicate entries by scanning the database for purely unique entries
### and collapsing those that are duplicates

App_Data <- Combined_Taxonomy %>%
  count(
    Dataset, Citation, Author, Family, Gene, Mutation, country, Source_Type, Host,
    name = "Unique_Entries"
  )


world <- ne_countries(scale = "medium", returnclass = "sf") %>%
  filter(name != "Antarctica")  

#-------------------------------------------------------------------------------
# Theme 
#-------------------------------------------------------------------------------

gene_colors <- c("#4b3b78", "#7b6bb0", "#5b8fc9", "#4fa3a5","#7fb77e",
                 "#e0b04f", "#d9825b", "#C8607a", "#9a9a9a" )

submission_template <- data.frame(
  'Family/species' = "Aspergillaceae", Gene = "cyp51a", Mutation = "g448S", Location = "USA",
  Accession = "MN0000", Sequence = "", Citation = "", AUthor = "", check.names = FALSE
)
required_cols <- c("Families/species", "Gene", "Mutation", "Location", "Source_Type","Accession")

app_theme <- bs_theme( 
  base_font = font_google("Inter"), primary = "lavender" , bg = "#FAFAF7", fg = "#1f2328" ) |>
  bs_add_rules(".value-box-value {font-size: 0.8rem; } .value-box-title { font-size: 0.5rem;}")

#-------------------------------------------------------------------------------
# UI - necessary for app deployment
#-------------------------------------------------------------------------------

ui <- page_navbar ( 
  header = useShinyjs(),
  title = tags$span(
    tags$img(src = "FELLogo.png", height = "20px", alt = "FEL logo", "  S.P.O.R.E The Global
           Fungal Antimicrobial Resistance Tracker")),
  window_title = "SPORE - Global Fungal Antimicrobial Resistance Tracker",
  theme = app_theme,
  nav_panel( "Home",
             tags$h2("Welcome to S.P.O.R.E", class = "text-center"),
             tags$p("This app allows you to explore current findings on fungal antimicrobial
                    resistance around the world, and even submit your own findings.", 
                    class = "text-center"),
             layout_columns(
               height = "100px",
               value_box(title = "Records", value = textOutput("n_records"),),
               value_box(title = "Countries", value = textOutput("n_countries"),),
               value_box(title = "Families", value = textOutput("n_families"))),
             maplibreOutput("home_globe", height = "650px")),
  nav_spacer(),
  nav_menu("Explore", 
           nav_panel("Mutations",
                     navset_tab(
                       nav_panel("Table", DTOutput("mutation_table")),
                       nav_panel("Top 20 mutations", plotOutput("mutation_bar", height = "600px"))
                     )
           ),
           nav_panel("Genes",
                     navset_tab(
                       nav_panel("Overview",
                                 selectInput("gene_picker","Choose a gene",
                                             choices = NULL),
                                 
                                 layout_columns(height = "200px",
                                                value_box(title = "Records", value = textOutput("gene_records")),
                                                value_box(title = "Families", value = textOutput("gene_families")),
                                                value_box(title = "Countries", value = textOutput("gene_countries")),
                                                value_box(title = "Mutations", value = textOutput("gene_mutations")))),
                       nav_panel("Table", DTOutput("gene_table")),
                       nav_panel("Records", plotOutput("gene_bar", height = "680px"))
                     )),
           nav_panel("Families",
                     navset_tab(
                       nav_panel("Alluvial Chart",
                                 selectInput("family_picker", "choose a family", choices = NULL), plotOutput("family_alluvial", height = "600px")),
                       layout_columns(height = "100px", value_box(title = "Records", value = textOutput("fam_records")),
                                      value_box(title = "Countries", value = textOutput("fam_countries")),
                                      value_box(title = "Genes", value = textOutput("fam_genes")),
                                      value_box(title = "Mutations", value = textOutput("fam_mutations"))),
                       nav_panel("Table",
                                 DTOutput("family_table")),
                       nav_panel("Records",plotOutput("family_bar", height = "600px")))),
           nav_panel("Add Data",
                     h3("Submit your data"),
                     p("Download our template, fill it out with one row per mutation record, then upload it below. We will review submissions and add them ASAP"),
                     downloadButton("download_template", "Download_template"),
                     fileInput("upload_file", "Upload template(.csv)", accept = ".csv"),
                     DTOutput("upload_preview"),
                     actionButton("submit_upload", "Submit", class = "btn-primary")), 
           nav_panel("Distribution", radioButtons("data_selection",label = NULL,
                                                  choices = c( Families = "Family",
                                                               Genes = "Gene",
                                                               Mutations = "Mutation"),
                                                  inline = TRUE),
                     girafeOutput("distribution_map"), uiOutput("lasso_check")),
           nav_panel("Methods",
                     h1("database"),
                     p("The S.P.O.R.E database was created in an attempt to track fungicide resistance on the global scale.
                       first a base database was created by using the FungAMR database and adding geodata to it that was extracted using provided accession numbers.
                       To do this a python script was written that would utilize each provided accession to search for associated geographic metadata.
                       This data was then returned into a new excel file, and any entry that did not have associated geodata was removed. Then
                       the database was fine tuned, seperating entries that had multiple genes or mutations per entry. Following this the authors read each entry
                       individually and removed entries where the geogeaphic data was misattributed either due to the accession being a reference, or the
                       study listing multiple samples from multiple countries without a concrete way to parse where each isolate came from. Reference accessions were identified 
                       by comparing the authors from the database entry to the authors from the cited papers and looking for a mismatch.
                       After this the database was expanded by the authors who manually added extra entries.
                       Finally a meta analysis was conducted to analyze patterns present in the data. This is where a unique entry was defined as a singular combination of 
                       mutation,Gene, and Family, which is stored on the website as records."),
      
           )), 
  
  nav_spacer(),
  nav_item(popover(
    tags$span(bs_icon("funnel"), "Filter data"),
    radioButtons(inputId = "source", label = NULL, choices = c(combined = "all", Agricultural = "Agricultural", Clinical = "Clinical"),
                 selected = "all", inline = TRUE),
    
    radioButtons( inputId = "database", label = NULL, choices = c(All = "all", FungAMR = "FAMR", FEL = "FEL"),
                  selected = "all", inline = TRUE))))






#-------------------------------------------------------------------------------
#server - necessary for app deployment
#-------------------------------------------------------------------------------

server <- function(input, output, session) {
  
  app_data<-reactive({
    out<-App_Data
    if(input$source != "all") out <- filter(out, Source_Type == input$source)
    if(input$database != "all") out <- filter(out, Dataset == input$database)
    out
  })
  
  output$n_records <- renderText(nrow(app_data()))
  
  output$n_countries <- renderText(n_distinct(app_data()$country, na.rm = TRUE))
  
  output$n_families <- renderText(n_distinct(app_data()$Family, na.rm = TRUE))
  
  country_summary <- reactive({
    app_data() %>%
      filter(!is.na(country)) %>%
      group_by(country) %>%
      summarise(
        families = n_distinct(Family),
        family_list = paste(sort(unique(Family)), collapse = ", "),
        popup = paste0("<b>", first(country),"</b><br>",
                       n_distinct(Family), " families: ",
                       paste(sort(unique(Family)), collapse = ", "))
      )
  })
  
  data_selection <- reactive({
    col <- input$data_selection
    label <- c(Family = "families", Gene = "genes", Mutation = "mutations") [col]
    summary <- app_data() %>%
      filter(!is.na(country)) %>%
      group_by(country) %>%
      summarise(
        count = n_distinct(.data[[col]]),
        popup = paste0("<br", first(country), "</b><b>",
                       n_distinct(.data[[col]])," ", label, ": ",
                       paste(sort(unique(.data[[col]])), collapse = ", "))
      )
    world %>% left_join(summary, by = c("name" = "country"))
    
  })
  
  map_data <- reactive({ 
    world %>%
      left_join(country_summary(), by = c("name" = "country"))
  })
  
  output$home_globe <- renderMaplibre({ 
    maplibre(
      style = carto_style("positron"),
      center = c(-30,20),
      zoom = 1.5 ) |>
      add_fill_layer (
        id = "country",
        source = map_data(),
        fill_color = interpolate(
          column = "families",
          values = c(1, max(map_data()$families, na.rm = TRUE)),
          stops = c("#E9E4F5","#4b3b78"), na_color = "#ececec"), fill_opacity = 0.5, fill_outline_color = "white",
        popup = "popup") |>
      add_legend("Fungal families", values = c(1, max(map_data()$families, na.rm = TRUE)),
                 colors = c("#E9E4F5","#4b3b78"), position = "bottom-left")|>
      htmlwidgets::onRender("
      function(el) {
        var map = el.map;
        if(window.matchMedia('(prefers-reduced-motion:reduce)').matches) return;
           var busy = false;
           function spin() { if (busy) return;
             var c = map.getCenter();
               c.lng -= 2;
               map.easeTo({center: c, duration:1000, easing: function(t) {return t;}});
      }
    map.on('mousedown', function() {busy = true;});
    map.on('touchstart', function() {busy = true;}) ;
    map.on('mouseup', function() {busy = false; spin();});
    map.on('touchend', function() {busy = false; spin();});
    setInterval(function(){if (!busy && !map.isMoving()) spin()}, 100);
    map.on('moveend', spin);
    spin();
    
    }")})
  
  output$distribution_map <- renderGirafe({ 
    p <-  ggplot(data_selection()) + geom_sf_interactive(aes(fill = count, tooltip = popup, data_id = name)) + coord_sf(crs = "+proj=eqearth") +
      theme(panel.background = element_rect(fill = "#eef3f8")) + scale_fill_gradient(low = "#E9E4F5", high = "#4b3b78", na.value = "#ececec", breaks = function(x) unique(round(pretty(x))))+ labs(fill = "Count")
    girafe(ggobj = p, width_svg = 11, height_svg = 5)
  })
  
  output$lasso_check <- renderUI({ sel <- input$distribution_map_selected 
  req(sel)
  label <- c(Family = "Families", Gene = "Genes", Mutation = "Mutations")[input$data_selection]
  fams <- app_data() %>%
    filter(country %in% sel) %>%
    pull(.data[[input$data_selection]]) %>%
    unique() %>%
    sort()
  tagList(
    p(strong("Countries selected: "), paste(sel, collapse = ", ")),
    p(strong("Present in selection: "), paste(fams, collapse = ", ")),
    p("Data Selected: ", input$data_selection)
  )})
  
  
  observe({ 
    updateSelectInput(session, "family_picker", choices = sort (unique(app_data()$Family)))
  })
  family_data <- reactive({
    req(input$family_picker)
    app_data() %>% filter (Family == input$family_picker)
  })  
  output$fam_records <- renderText(nrow(family_data()))
  output$fam_countries <- renderText(n_distinct(family_data()$country, na.rm = TRUE))
  output$fam_genes <- renderText(n_distinct(family_data()$Gene, na.rm = TRUE))
  output$fam_mutations <- renderText(n_distinct(family_data()$Mutation, na.rm = TRUE))
  output$family_alluvial <- renderPlot({
    req(input$family_picker)
    plot_data <- app_data() %>%
      filter(Family == input$family_picker, !is.na(country)) %>%
      count(Family, country, Gene, Mutation)
    
    output$family_table <- renderDT({
      app_data()%>%
        group_by(Family) %>%
        summarise(Records = n(),
                  Countries = paste(sort(unique(na.omit(country))), collapse = ", "),
                  Genes = paste(sort(unique(na.omit(Gene))), collapse = ", ")) %>%
        datatable( rownames = FALSE)
    })
    
    output$family_bar <- renderPlot({
      app_data()%>%
        count(Family) %>%
        ggplot(aes(x = n, y = reorder(Family,n))) +
        geom_col(fill = "#4b3b78") +
        labs(x = "Records", y = NULL)
    })
    
    ggplot(plot_data, aes(axis1 = Family, axis2 = country, axis3 = Gene, axis4 = Mutation, y = n)) +
      geom_alluvium(aes(fill = Gene)) +
      geom_stratum() +
      geom_text(stat = "stratum", aes(label = after_stat(stratum)), size =3) +
      scale_x_discrete(limits =c("Family","Country","Gene","Mutation"))+ theme(panel.background =  element_blank(),
                                                                               axis.text.y = element_blank(),
                                                                               axis.ticks.y = element_blank()) +
      scale_fill_manual(values = gene_colors) + labs(y=NULL) 
  })
  
  
  observe({ updateSelectInput(session, "gene_picker",
                              choices = sort(unique(app_data()$Gene)))})
  gene_data <- reactive({
    req(input$gene_picker)
    app_data() %>% filter(Gene == input$gene_picker)
  })
  output$gene_records <- renderText(nrow(gene_data()))
  output$gene_families <- renderText(n_distinct(gene_data()$Family, na.rm =TRUE))
  output$gene_countries <- renderText(n_distinct(gene_data()$country, na.rm = TRUE))
  output$gene_mutations <- renderText(n_distinct(gene_data()$Mutation, na.rm = TRUE))
  
  output$gene_table <-renderDT({
    app_data() %>%
      group_by(Gene) %>%
      summarise(Records = n(),
                Families = paste(sort(unique(na.omit(Family))), collapse = ", "),
                Countries = paste(sort(unique(na.omit(country))), collapse = ", ")) %>%
      datatable(rownames = FALSE)
  })
  output$gene_bar <- renderPlot({
    app_data() %>%
      count(Gene) %>%
      ggplot(aes(x = n, y = reorder(Gene, n))) +
      geom_col(fill = "#4b3b78") +
      labs(x = "Records", y = NULL)
  })
  
  output$mutation_table <- renderDT({
    app_data() %>%
      group_by(Mutation) %>%
      summarise(Records = n(),
                Genes = paste(sort(unique(na.omit(Gene))), collapse = ", "),
                Families = paste(sort(unique(na.omit(Family))), collapse = ", "),
                Countries = paste(sort(unique(na.omit(country))), collapse = ", ")) %>%
      arrange(desc(Records)) %>%
      datatable(rownames = FALSE)
  })
  output$mutation_bar <- renderPlot({
    app_data() %>%
      count(Mutation) %>%
      slice_max(n, n = 20, with_ties = FALSE) %>%
      ggplot(aes(x=n, y = reorder(Mutation, n))) +
      geom_col(fill = "#4b3b78") +
      labs(x = "records", y = NULL)
  })
  
  output$download_template <- downloadHandler(
    filename = "SPORE_submittion_template.csv",
    content = function(file) write.csv(submission_template, file, row.names = FALSE)
  )
  upload_data <- reactive({
    req(input$upload_file)
    df <- tryCatch(read.csv(input$upload_file$datapath, check.names = FALSE,
                            colClasses = "character"),
                   error = function(e) NULL)
    validate(need(!is.null(df), "Please upload data as a .csv"))
    validate(need(nrow(df) > 0, "Please fill in data" ))
    missing <- setdiff(required_cols, names(df))
    validate(need(length(missing) == 0, paste("missing required data:", paste(missing, collapse = ", "))))
    df[] <- lapply(df, trimws)
    blank_rows <- which(rowSums(df[required_cols] == "" |is.na(df[required_cols])) > 0)
    validate(need(length(blank_rows) == 0,
                  paste("required cells are empty in:",
                        paste(blank_rows + 1, collapse = ", "))))
    validate(need(!any(df$Accession == "e.g MN00000"),
                  "Please replace ot deleate example row"))
    validate(need(!any(grepl("^[=@]",unlist(df)))),
             "cells can't start with = or @, try again")
    df
  })
  output$upload_preview <- renderDT(datatable(upload_data(), rownames = FALSE))
  
  
  observeEvent(input$submit_upload, {
    req(upload_data())
    tmp <- tempfile(fileext = ".csv")
    write.csv(upload_data(), tmp, row.names = FALSE)
    path <- paste0("submissions/", format(Sys.time(), "%Y-%m-%d_%H%M%S"),
                   "_",sample(1000:9999,1), ".csv")
    tryCatch({ gh::gh("PUT /repos/{owner}/{repo}/contents/{path}",
                      owner = "kcrovetti1", repo = "spore-submissions", path=path,
                      message = paste("NEW submission", path),
                      content = jsonlite::base64_enc(readBin(tmp, "raw",
                                                             file.size(tmp))),
                      .token = Sys.getenv("SPORE_GITHUB_PAT"))
      showNotification("Thanks! we will review your submission as soon as we can!", type = "message")},
      error = function(e) {
        message("upload failed: ", conditionMessage(e))
        showNotification("Submission failed, try again later",
                         type = "error")
      })
  })
  
}
#-------------------------------------------------------------------------------
#deployment
#-------------------------------------------------------------------------------

shinyApp(ui,server)



