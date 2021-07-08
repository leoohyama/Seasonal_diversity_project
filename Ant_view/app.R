library(shiny)
library(ggplot2)
df<-readRDS("data_ants.rds")
species_traits<-read.csv("species_traits_florida.csv")



ui <- fluidPage(
    selectInput("Speciesone", "Species 1", unique(df$Species)),
    selectInput("Speciestwo", "Species 2", unique(df$Species)),
    plotOutput("plot_ant"),
    plotOutput("webplot"),
    plotOutput("plot_ant1"),
    plotOutput("webplot1")
)

server <- function(input, output, server) {
    
    dat<-reactive({
        dplyr::filter(df, Species == input$Speciesone)
    })
    sp1<-reactive({species_traits %>% dplyr::filter(species == input$Speciesone) %>%
        dplyr::select(species,WL, HW, EL, ML) %>%
        pivot_longer(cols = c(WL, HW, EL, ML))})
    
    datone<-reactive({
        dplyr::filter(df, Species == input$Speciestwo)
    })
    sp1one<-reactive({species_traits %>% dplyr::filter(species == input$Speciestwo) %>%
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
                  
                  strip.text = element_text(family = "Helvetica", face= "bold",
                                            size =12))
    },
    height = 400,width = 600)
    
    output$webplot<-renderPlot({
        
        ggplot(data = sp1(), aes(x = name, y = value)) +
            geom_bar(stat = "identity") + coord_polar()
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
                  
                  strip.text = element_text(family = "Helvetica", face= "bold",
                                            size =12))
    },
    height = 400,width = 600)
    
    output$webplot1<-renderPlot({
        
        ggplot(data = sp1one(), aes(x = name, y = value)) +
            geom_bar(stat = "identity") + coord_polar()
    },
    height = 400,width = 600)
}

shinyApp(ui, server)




