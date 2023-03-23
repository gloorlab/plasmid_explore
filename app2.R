library(shiny)

# Define UI for application that draws a histogram
ui <- fluidPage(

    
    # Application title
    titlePanel("data stripchart"),
    
    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            # set a directory name
            # input$upload will have the default of bakta
            
            sliderInput("pch", # choose a symbol, input$pch
                        "symbol choice:",
                        min = 0,
                        max = 20,
                        value = 19),
            selectInput("sample", "choose sample",
                        c("1-T1"="1-T1",
                          "2-T4"='2-T4',
                          "3-T4",'3-T4',
                          "4-T3", "4-T3",
                          "5-T2",'5-T2'
                            ))
        ),
        
        # Show a plot of the generated distribution
        mainPanel(
            plotOutput("distPlot")
            #h4("output$dir"),
            #verbatimTextOutput("empty")
        )
    )
)


# Define server logic required to draw a stripchart
server <- function(input, output, session) {
    
    output$distPlot <- renderPlot({
        
        # the dir can be in a dropdown box
        # this is called every time the gb is called to open the file
        up.data <- reactive({
              dir <- paste('~/Documents/thesis/data/', input$sample,'/', sep='')
        })
        
        x <- list.files(up.data(), pattern='[0-9].tsv') # tsv file only
        # move up.files to a simple variable for convenience
        
        annot <- matrix(data=NA, nrow=length(x), ncol=5)
        len <- vector()
        keggFound <- vector()
        
        for(i in 1:length(x)){
            gb <- read.table(paste(up.data(),x[i], sep=''), comment.char="#", sep='\t', quote='')
            annot[i,1] <- dim(gb)[1] # total cds
            annot[i,2] <- length(grep('PFAM', gb[,9])) # PFAM annotated ORFS in column 9
            kolines <- grep('KEGG', gb[,9])
            for(k in kolines){ keggFound <- c(keggFound, gsub(".*KEGG:(K\\d+).+", "\\1", gb[k,9]) )}
            annot[i,4] <- length(kolines) # KEGG annotated ORFS in column 9
            annot[i,3] <- annot[i,2]/annot[i,1] # proportion n.annot / n.rows
            annot[i,5] <- annot[i,4]/annot[i,1] # proportion n.annot / n.rows
            len <- c(len, gb[,4] - gb[,3])
        }
        
        K0_beta <- c("K18093","K09476","K09475","K08720","K18133","K12552","K05366","K03693","K12555","K00687","K12556","K05515","K12553","K03587","K18129","K18130","K18131","K18135","K18136","K18137","K18143","K18144")
        K0_napt <- c("K14579","K14580","K14578","K14581","K14582","K14583","K14584","K14585","K00152","K18242","K18243","K14578","K14581")
        K0_pc <- c("K18074","K18075","K18077","K18076","K18068","K18069","K18076","K04102","K18251","K180252","K180253 ","k180254","K18256")
        overlapBeta <- intersect(keggFound, K0_beta)
        overlapNapt <- intersect(keggFound, K0_napt)
        overlapPc <- intersect(keggFound, K0_pc)
        olB.name <- paste("KO_beta:", overlapBeta ,sep=',') 
        olN.name <- paste("KO_napt", overlapNapt ,sep=',') 
        olP.name <- paste("KO_napt", overlapPc ,sep=',') 
        glab <- c('PFAM=1', 'KEGG=2', 'K0_beta')
        glab <- c('PFAM=1', 'KEGG=2', 'K0_napt')
        glab <- c('PFAM=1', 'KEGG=2', 'K0_pc')
        
        # plot the stripchart with the column and symbol of choice
        # in rgb, the last value is the transparency from 0-1
        # stripchart(list(annot1[,3], annot2[,3]), rest of parameters )
        stripchart(list(annot[,3], annot[,5], 0), method='jitter', group.names=glab, vertical=T, pch=input$pch, col=rgb(1,0,0,0.2))
        text(3,max(annot[,3]), labels=olB.name)
        text(3,max(annot[,4]), labels=olN.name)
        text(3,max(annot[,5]), labels=olP.name)
        
      })
    
}

# Run the application 
shinyApp(ui = ui, server = server)