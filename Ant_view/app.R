library(shiny)
library(tidyverse)
df<-readRDS("data_ants.rds")
species_traits<-read.csv("species_traits_florida.csv")
library(patchwork)



ui <- fluidPage(
    
    
    titlePanel("Compare traits and seasonal dynamics of the ants of Central Florida!"),
    tableOutput("data"),
    fluidRow(
        column(4,selectInput("Speciesone", "Species 1", unique(df$Species))),
        column(4,selectInput("Speciestwo", "Species 2", unique(df$Species)))
        
    ),
    plotOutput("plot_ant", width = "100%", click = "plot_click"),
    plotOutput("plot_ant1", width = "100%", click = "plot_click"),
    plotOutput("plot_ant2", width = "100%"),
    

    
)

server <- function(input, output, server) {

    
    dat<-reactive({
        dplyr::filter(df, Species == input$Speciesone)
    })
    datone<-reactive({
        dplyr::filter(df, Species == input$Speciestwo)
    })
    sp1<-reactive({species_traits %>% dplyr::filter(species %in% c(input$Speciesone, input$Speciestwo)) %>%
        dplyr::select(species,WL, HW, EL, ML) %>%
        pivot_longer(cols = c(WL, HW, EL, ML))})
    
    
    
    
    output$plot_ant<-renderPlot({
        ggplot(dat()) + 
            geom_rect(xmin=17652,
                      xmax=17806, ymin=-Inf, ymax=Inf, fill="lightblue", alpha=1) +
            geom_rect(xmin=18017,
                      xmax=Inf, ymin=-Inf, ymax=Inf, fill="lightblue", alpha=1) +
            geom_line(aes(x = time_code, y = (total_f), color = Species, group = Species), size = 2, color = "darkgreen") +
            geom_point(aes(x = time_code, y = (total_f)),
                       pch =21, size = 4, color = "black", fill = "darkgreen") +
            labs(y = "Abundance", x = "Time", title = input$Speciesone) + theme_dark() +
            theme(axis.text.x = element_text(angle = 45, face = "bold", size = 12, vjust =  0.5),
                  legend.position = "none" ,
                  axis.ticks = element_blank(),
                  plot.margin = unit(c(0,0,0,0), "cm"),
                  strip.text = element_text(family = "Helvetica", face= "bold",
                                            size =12))
        
     
    },
    height = 400,width = 600)
    
    
    output$plot_ant1<-renderPlot({
        ggplot(datone()) + 
            geom_rect(xmin=17652,
                      xmax=17806, ymin=-Inf, ymax=Inf, fill="lightblue", alpha=1) +
            geom_rect(xmin=18017,
                      xmax=Inf, ymin=-Inf, ymax=Inf, fill="lightblue", alpha=1) +
            geom_line(aes(x = time_code, y = (total_f), color = Species, group = Species), size = 2, color = "darkgreen") +
            geom_point(aes(x = time_code, y = (total_f)),
                       pch =21, size = 4, color = "black", fill = "darkgreen") +
            labs(y = "Abundance", x = "Time", title = input$Speciestwo) + theme_dark() +
            theme(axis.text.x = element_text(angle = 45, face = "bold", size = 12, vjust =  0.5),
                  legend.position = "none" ,
                  axis.ticks = element_blank(),
                  plot.margin = unit(c(0,0,0,0), "cm"),
                  strip.text = element_text(family = "Helvetica", face= "bold",
                                            size =12))
    },height = 400,width = 600)
    
    output$plot_ant2<-renderPlot({
        ggplot(data = sp1(), aes(x = name, y = value, fill = species)) +
            geom_bar(stat = "identity") + coord_polar() +
            labs(title = "Trait values", y = "Measurement in milimeters", x = NA) + 
            theme_dark() +
            theme(plot.margin = unit(c(0,0,0,0), "cm"))
    },height = 400,width = 600)
    
    output$data <- renderTable({
        nearPoints(dat(), input$plot_click, xvar = "time_code", yvar = "total_f")
    })
    
}

shinyApp(ui, server)




