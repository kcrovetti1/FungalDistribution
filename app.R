#README:

#Hello, thank you for checking out the code behind our app. I started the app for comparing
#Some data for one of my project and as it great so did the functions of the app
#I have gone back and cleaned up the app to make the code as understandable as possible
#So that anyone can understand and make edits. The code was written by myself through trial
# And error supported with a lot of googling and reading stack overflow. Because I only
#Started learning R 7 months ago it is fairly rudimatary with several repeating scrtips
#Feel free to recomend ways to streamline the app functions. Also before publishing the app 
#ChatGPT was used to stresstest and debug areas of the code but its involvement was limited to this.
#I hope you enjoy the app and learn something from the data!

#-------------------------------------------------------------------------------
#Load dependencies 
#-------------------------------------------------------------------------------
library(shiny)
library(bslib)
library(dplyr)
library(ggplot2)
library(rnaturalearth)
library(sf)
library(stringr)
library(tidyr)
library(leaflet)
library(purrr)
library(scales)
library(DT)
library(ggalluvial)
library(readxl)
library(viridis)
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

FungAMR_Combined <- bind_rows(
  FungAMR_Clinical %>% mutate(origin = "Clinical"),
  FungAMR_Agricultural %>% mutate(origin = "Agricultural"))
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
# UI
#-------------------------------------------------------------------------------
ui <- page_sidebar(
  title = tags$span(
    tags$img(
      src = "FELLogo.png",
      height = "40px",
      style = "margin-right: 15px; vertical-align: middle;"
    ),
    tags$span(
      "FungAMR & FEL Database Comparison",
      style = "vertical-align: middle; font-size: 24px;"
    )
  ),
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
        "Clinical Data"           = "clinical",
        "Agricultural Data"       = "agricultural",
        "Resistance Genes"        = "genes",
        "Mutations"               = "mutations",
        "Geographic Distribution" = "geographic",
        "Data by Family"          = "family",
        "Compare Databases"       = "compare"
      ),
      selected = "home"
    )), uiOutput("page_content"))
#-------------------------------------------------------------------------------
# Server
#-------------------------------------------------------------------------------
server <- function(input, output, session) {
  
  make_count_table <- function(df, col, label, page_length = 10) {
    tryCatch({
      result <- df %>%
        filter(!is.na(.data[[col]])) %>%
        count(!!sym(label) := .data[[col]], sort = TRUE)
      
      datatable(
        result,
        options = list(pageLength = page_length)
      )
    }, error = function(e) {
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
        families_display = case_when(
          is.na(families) ~ "No families recorded",
          families == "" ~ "No families recorded",
          TRUE ~ families
        ),
        
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
      addControl(
        html = paste0(
          "<div style='background:white;padding:8px 12px;",
          "border-radius:5px;font-size:13px;'>",
          
          "<div style='margin-bottom:5px;'>",
          "<span style='display:inline-block;width:14px;height:14px;",
          "background:#08519C;margin-right:6px;vertical-align:middle;'></span>",
          "Data present",
          "</div>",
          
          "<div>",
          "<span style='display:inline-block;width:14px;height:14px;",
          "background:#BDBDBD;margin-right:6px;vertical-align:middle;'></span>",
          "No data",
          "</div>",
          
          "</div>"
        ),
        position = "bottomright"
      )
  }
  

  
  output$page_content <- renderUI({
    switch(
      input$page_selection,
      

      
      "home" = div(
        h2("Welcome to the Database Comparison Tool"),
        p("This Shiny app allows you to explore and compare FungAMR and FEL databases."),
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
        )
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
            "Genes Comparison",
            p("FungAMR Clinical Genes:"),
            DTOutput("fungamr_clinical_genes"),
            br(),
            p("FEL Clinical Genes:"),
            DTOutput("fel_clinical_genes")
          ),
          
          tabPanel(
            "Mutations Comparison",
            p("FungAMR Clinical Mutations:"),
            DTOutput("fungamr_clinical_mutations"),
            br(),
            p("FEL Clinical Mutations:"),
            DTOutput("fel_clinical_mutations")
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
            "Genes Comparison",
            p("FungAMR Agricultural Genes:"),
            DTOutput("fungamr_agricultural_genes"),
            br(),
            p("FEL Agricultural Genes:"),
            DTOutput("fel_agricultural_genes")
          ),
          
          tabPanel(
            "Mutations Comparison",
            p("FungAMR Agricultural Mutations:"),
            DTOutput("fungamr_agricultural_mutations"),
            br(),
            p("FEL Agricultural Mutations:"),
            DTOutput("fel_agricultural_mutations")
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
            "Comparison",
            
            h4("FungAMR vs FEL Resistance Genes"),
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
            "Comparison",
            
            h4("FungAMR vs FEL Resistance Mutations"),
            p(
              "Comparison of resistance-associated mutations identified in each database."
            ),
            
            DTOutput("resistance_mutation_comparison")
          )
        )
      ),
      
     
      
      "geographic" = div(
        h2("Geographic Distribution(Please Give Time To Load"),
        
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
            
            br(),
            
            h4("Agricultural Records"),
            
            p(
              "FungAMR Agricultural: ",
              nrow(FungAMR_Agricultural)
            ),
            
            p(
              "FEL Agricultural: ",
              nrow(FEL_Agricultural)
            )
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
      )
    )
  })
  

  output$fungamr_clinical_table <- renderDT(
    make_count_table(
      FungAMR_Clinical,
      "Family",
      "Family"
    )
  )
  
  output$fel_clinical_table <- renderDT(
    make_count_table(
      FEL_Clinical,
      "Family/species",
      "Family"
    )
  )
  
  output$fungamr_clinical_genes <- renderDT(
    make_count_table(
      FungAMR_Clinical,
      "gene.or.protein.name",
      "Gene"
    )
  )
  
  output$fel_clinical_genes <- renderDT(
    make_count_table(
      FEL_Clinical,
      "Gene",
      "Gene"
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
  
  output$fel_agricultural_genes <- renderDT(
    make_count_table(
      FEL_Agricultural,
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
    
    comparison <- full_join(
      fungamr_genes,
      fel_genes,
      by = "Gene"
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
        Total = FungAMR_Count + FEL_Count,
        Difference = FungAMR_Count - FEL_Count,
        Database = case_when(
          FungAMR_Count > 0 &
            FEL_Count > 0 ~ "Both",
          FungAMR_Count > 0 ~ "FungAMR Only",
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
    
    comparison <- full_join(
      fungamr_mutations,
      fel_mutations,
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
        Total = FungAMR_Count + FEL_Count,
        Difference = FungAMR_Count - FEL_Count,
        Database = case_when(
          FungAMR_Count > 0 &
            FEL_Count > 0 ~ "Both",
          FungAMR_Count > 0 ~ "FungAMR Only",
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
    
    all_families <- c(
      fungamr_families,
      fel_families
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
    
    comparison <- full_join(
      fungamr_genes,
      fel_genes,
      by = c("Gene")
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
    
    comparison <- full_join(
      fungamr_mutations,
      fel_mutations,
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

#-------------------------------------------------------------------------------
# Run
#-------------------------------------------------------------------------------
shinyApp(
  ui,
  server
)