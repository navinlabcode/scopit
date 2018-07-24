library(shiny)
library(ggplot2)
library(pmultinom)
library(stringr)

source("backend.R")

# Display the error messages that I give, not a generic error message
options(shiny.sanitize.errors = FALSE)

aislyn.palette <- c("tomato2","tan2","darkseagreen3","cadetblue","royalblue4","grey80","grey55")

make.minus <- function(i, input, session, clone.to.remove, current.clones)
{
    print(sprintf("Listening for %d", i))
    observeEvent(input[[sprintf("minus%d", i)]], {
      print(sprintf("Removing clone %d", i))
      # The clone frequencies draws input from boxes which are defined in clonefreq.boxes
      # Both those boxes and numeric.clonefreqs will respond to a change in the number of clones
      # When this reaction happens, it has to be recorded which has to be removed
      clone.to.remove(i)
      # Decrement the number of clones
      old.value <- current.clones()
      print(sprintf("Currently there's %d clones", old.value))
      # Trying to put in a hack so that new clones will always have the last inputted frequency
      # If what you're minusing is the bottom clone, and its not the only clone,
      # set its frequency to the frequency of the previous clone before deleting
      if (i == old.value & i > 1)
      {
        prev.freq <- input[[sprintf("clonefreq%d", i-1)]]
        this.id <- sprintf("clonefreq%d", i)
        # I have no idea why I have to make this reactive
        # BUt if I just put updateNumericInput here, it won't happen
        observeEvent(current.clones(), {
          print(sprintf('Resetting clone "%s" frequency to %f', this.id, prev.freq))
          updateNumericInput(session, this.id, value=prev.freq)
        }, once=TRUE, ignoreInit=TRUE)
      }
      current.clones(old.value-1)
      print(sprintf("Now there's %d", current.clones()))
    }, once=TRUE, ignoreInit=TRUE)
}
ui <- fluidPage(
  titlePanel("Number of cells to sequence, v1.0.0"),
  p("by Alexander Davis"),
  tabsetPanel(tabPanel(title="Prospective",
  sidebarLayout(sidebarPanel(id="sidebar",
                             numericInput("cutoff", "Number of cells of each type which must be sequenced", 3),
                             numericInput("power", "Required probability of sequencing this many cells from each type", .99,
                                          min=0, max=1, step=.01),
                             uiOutput("clonefreqs"),
                             actionButton("plus_button", "+")),
                mainPanel(plotOutput("path_plot"))
  )
  )),
                theme="theme.css"
)

server <- function(input, output, session)
{
  current.clones <- reactiveVal(1)
  
  # Keep track of which clone to remove
  clone.to.remove <- reactiveVal(NA)
  
  # Make an observer on the first minus button
  make.minus(1, input, session, clone.to.remove, current.clones)

  # Respond to the plus button
  observeEvent(input$plus_button, {
    old.value <- current.clones()
    new.value <- old.value + 1
    current.clones(new.value)
    
    make.minus(new.value, input, session, clone.to.remove, current.clones)
  })
  
  numeric.clonefreqs <- reactive({
    these.clonefreqs <- sapply(1:current.clones(), function(i) input[[sprintf("clonefreq%d",i)]])
    these.clonecounts <- sapply(1:current.clones(), function(i) input[[sprintf("clonecount%d",i)]])
    if (any(sapply(these.clonefreqs, is.null))) {
      NULL
    } else {
      rep.int(these.clonefreqs, these.clonecounts)
    }
    })
  numeric.path <- reactive({
    if (any(numeric.clonefreqs() > 1)) stop("Subpopulation frequencies must be less than 1")
    if (sum(numeric.clonefreqs()) > 1) stop("Total subpopulation frequency must be less than 1")
    if (any(numeric.clonefreqs() < 0)) stop("Subpopulation frequencies cannot be negative")
    if (current.clones() < 1) stop("Need at least one subpopulation")
    numeric.allfreqs <- c(numeric.clonefreqs(), 1 - sum(numeric.clonefreqs()))
    # You have to subtract 1 from the cutoff to get a maximum excluded
    all.cutoffs <- c(rep.int(input$cutoff-1, length(numeric.clonefreqs())), -1)
    #numeric.path <- upto(input$power, 3, all.cutoffs, numeric.allfreqs)
    required.samplesize <- invert.pmultinom.stoppable(lower=all.cutoffs, probs=numeric.allfreqs, target.prob=input$power, method="exact", maxval=10001)
    if (!is.na(required.samplesize)) {
      numeric.path <- pmultinom(lower=all.cutoffs, size=0:required.samplesize, probs=numeric.allfreqs, method="exact")
    } else {
      numeric.path <- pmultinom(lower=all.cutoffs, size=0:10001, probs=numeric.allfreqs, method="exact")
    }
  })


  recnum <- reactive({
    # Recommended number of cells to sequence
    recnum <- which(numeric.path() >= input$power)[1] - 1
  })

  output$path_plot <- renderPlot({
    # Get the probability values that will be plotted
    this.numeric.path <- numeric.path()

    # Pair these probability values with their associated sample size
    plot.data <- data.frame(x=0:(length(this.numeric.path)-1), y=this.numeric.path)

    # Make the plot
    if (!is.na(recnum())) {
      this.title <- ggtitle(sprintf("%d cells required for %.2f probability of success", recnum(), input$power))
    } else {
      this.title <- ggtitle(sprintf("%d insufficient for %.2f probability of success", length(this.numeric.path)+1, input$power))
    }
    ggplot(data=plot.data, mapping=aes(x=x, y=y)) +
      geom_line(size=rel(1.3)) +
      geom_hline(yintercept=input$power, colour="tan2", size=rel(1.3)) +
      geom_vline(xintercept=recnum(), colour="tan2", size=rel(1.3)) +
      xlab("Number of cells") +
      ylab(sprintf("Probability of ≥%d from each subpopulation", input$cutoff)) +
      # ylim must be set very slightly below 0, so that 0 values with small
      # downward error will display
      ylim(-0.000001,1.000001) +
      theme(axis.title=element_text(size=rel(1.7)),
            axis.text=element_text(size=rel(1.5)),
            plot.title=element_text(hjust=.5, size=rel(1.7)),
            panel.background=element_blank(),
            axis.line=element_line(size=rel(1.3))) +
      this.title
  })

  clonefreq.boxes <- eventReactive(current.clones(), {
    # What is the "current" clone frequency--that is, the frequency of the bottom clone on the list?
    clone.freqs <- lapply(1:current.clones(), function(i) input[[sprintf("clonefreq%d",i)]])
    if (all(sapply(clone.freqs, is.null))) {
      current <- .01
    } else {
      current <- clone.freqs[[max(which(!sapply(clone.freqs, is.null)))]]
    }

    # Make a mapping from new numbering to old numbering
    this.clone.to.remove <- clone.to.remove()
    if (is.na(this.clone.to.remove)) {
      old <- 1:current.clones()
    } else {
      before.removal <- 1:(current.clones()+1)
      old <- before.removal[-this.clone.to.remove]
    }
    # Make the boxes that allow user input of clone frequencies
    if (current.clones() == 0) {
      clone.indices <- integer(0)
    } else {
      clone.indices <- 1:current.clones()
    }
    cancer.boxes <- lapply(clone.indices, function(i)
      flowLayout(
        numericInput(sprintf("clonefreq%d", i), ifelse(i==1, "Frequency of rarest type", "Frequency of additional types"),
                     ifelse(is.null(input[[sprintf("clonefreq%d",old[i])]]), current, input[[sprintf("clonefreq%d",old[i])]]), min=0.01, max=1, step=.01),
        numericInput(sprintf("clonecount%d", i), ifelse(i==1, "# of types with the lowest frequency", "# of types with this frequency"), 1),
        verticalLayout(actionButton(sprintf("minus%d", i), "-"))
      )
    )
    
    # Since the clone has already been removed, set clone.to.remove to NA
    clone.to.remove(NA)
    
    # Put these boxes into a vertical layout and return
    do.call(verticalLayout, cancer.boxes)
  })
  output$clonefreqs <- renderUI({
    clonefreq.boxes()
  })
}

shinyApp(ui=ui, server=server)
