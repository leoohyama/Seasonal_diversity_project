library(shiny)
library(ggplot2)

df<-readRDS("data_ants.rds")



ui <- fluidPage(
    selectInput("Species", "What's species do you want to see", unique(df$Species)),
    plotOutput("plot_ant")
)

server <- function(input, output, server) {
    
    dat<-reactive({
        dplyr::filter(df, Species == input$Species)
    })
    output$plot_ant<-renderPlot({
        ggplot(dat()) + 
            geom_rect(xmin=17652,
                      xmax=17806, ymin=-Inf, ymax=Inf, fill="lightblue", alpha=1) +
            geom_rect(xmin=18017,
                      xmax=Inf, ymin=-Inf, ymax=Inf, fill="lightblue", alpha=1) +
            geom_line(aes(x = time_code, y = (total_f), color = Species, group = Species), size = 2) +
            geom_point(aes(x = time_code, y = (total_f), fill = Species),
                       pch =21, size = 2, color = "black") +
            labs(y = "Abundance", x = "Time") + theme_dark() +
            theme(axis.text.x = element_text(angle = 45, face = "bold", size = 12, vjust =  0.5),
                  legend.position = "none" ,
                  axis.ticks = element_blank(),
                  
                  strip.text = element_text(family = "Helvetica", face= "bold",
                                            size =12))
    },
    height = 400,width = 600)
}

shinyApp(ui, server)



