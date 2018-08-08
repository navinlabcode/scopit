library(shiny)
library(ggplot2)
library(pmultinom)
library(stringr)
library(htmltools)
library(markdown)

source("backend.R")

# Display the error messages that I give, not a generic error message
options(shiny.sanitize.errors = FALSE)

make.minus <- function(i, input, session, clone.to.remove, current.clones, prefix, clonefreqname)
{
  print(sprintf("Minus button %s%d listening for %s%d", prefix, i, clonefreqname, i))
  observeEvent(input[[sprintf("%s%d", prefix, i)]], {
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
      prev.freq <- input[[sprintf("%s%d", clonefreqname, i-1)]]
      this.id <- sprintf("%sfreq%d", clonefreqname, i)
      # I have no idea why I have to make this reactive
      # But if I just put updateNumericInput here, it won't happen
      observeEvent(current.clones(), {
        print(sprintf('Resetting clone "%s" frequency to %f', this.id, prev.freq))
        updateNumericInput(session, this.id, value=prev.freq)
      }, once=TRUE, ignoreInit=TRUE)
    }
    current.clones(old.value-1)
    print(sprintf("Now there's %d", current.clones()))
  }, once=TRUE, ignoreInit=TRUE)
}

# Path calculation
# Making a function so I don't have to copy and paste for both panels
numeric.path.maker <- function(input, numeric.clonefreqs, current.clones, requested.point=NULL)
{
  reactive({
    if (any(numeric.clonefreqs() > 1)) stop("Type frequencies must be less than 1")
    if (sum(numeric.clonefreqs()) > 1) stop("Total type frequency must be less than 1")
    if (any(numeric.clonefreqs() < 0)) stop("Type frequencies cannot be negative")
    if (current.clones() < 1) stop("Need at least one type")
    numeric.allfreqs <- c(numeric.clonefreqs(), 1 - sum(numeric.clonefreqs()))
    # You have to subtract 1 from the cutoff to get a maximum excluded
    all.cutoffs <- c(rep.int(input$cutoff-1, length(numeric.clonefreqs())), -1)
    required.samplesize <- invert.pmultinom.stoppable(lower=all.cutoffs, probs=numeric.allfreqs, target.prob=input$power, method="exact", maxval=10001)
    if (!is.null(requested.point)) {
      right.limit <- min(max(required.samplesize, input[[requested.point]]), 10001)
    } else {
      right.limit <- min(required.samplesize, 10001)
    }
    if (!is.na(right.limit)) {
      numeric.path <- pmultinom(lower=all.cutoffs, size=0:right.limit, probs=numeric.allfreqs, method="exact")
    } else {
      numeric.path <- pmultinom(lower=all.cutoffs, size=0:10001, probs=numeric.allfreqs, method="exact")
    }
  })
}

# Function to produce a plot. This is really complicated and I don't want to
# copy and paste it for the second panel so I'm making a function.
path.plot.maker <- function(input, numeric.path, recnum, power, cutoff)
{
  # Get the probability values that will be plotted
  this.numeric.path <- numeric.path()
  
  # Pair these probability values with their associated sample size
  plot.data <- data.frame(x=0:(length(this.numeric.path)-1), y=this.numeric.path)
  
  # Make the plot
  if (!is.na(recnum())) {
    this.title <- ggtitle(sprintf("%d cells required for %.2f probability of success", recnum(), input[[power]]))
  } else {
    this.title <- ggtitle(sprintf("%d insufficient for %.2f probability of success", length(this.numeric.path)+1, input[[power]]))
  }
  ggplot(data=plot.data, mapping=aes(x=x, y=y)) +
    geom_line(size=rel(1.3)) +
    geom_hline(yintercept=input[[power]], colour="tan2", size=rel(1.3)) +
    geom_vline(xintercept=recnum(), colour="tan2", size=rel(1.3)) +
    xlab("Number of cells") +
    ylab(sprintf("Probability of ≥%d of each type", input[[cutoff]])) +
    # ylim must be set very slightly below 0, so that 0 values with small
    # downward error will display
    ylim(-0.000001,1.000001) +
    xlim(0,min(10001, max(plot.data$x))) + 
    theme(axis.title=element_text(size=rel(1.7)),
          axis.text=element_text(size=rel(1.5)),
          plot.title=element_text(hjust=.5, size=rel(1.7)),
          panel.background=element_blank(),
          axis.line=element_line(size=rel(1.3))) +
    this.title
}


# Basic settings, which will be part of both panels
basic.settings <- verticalLayout(numericInput("cutoff", "Number of cells of each type which must be sequenced", 3),
                                 numericInput("power", "Required probability of sequencing this many cells from each type", .95,
                                              min=0, max=1, step=.01))
# A clone of the basic settings, since Shiny doesn't support having the same element in two places
basic.settings.2 <- verticalLayout(numericInput("cutoff2", "Number of cells of each type which must be sequenced", 3),
                                   numericInput("power2", "Required probability of sequencing this many cells from each type", .99,
                                                min=0, max=1, step=.01))
# This will require code that keeps them in sync
ui <- fluidPage(
  titlePanel("SCOPIT, v1.0.0"),
  p("by Alexander Davis"),
  tabsetPanel(id="tab",
              tabPanel(title="Prospective",
                       sidebarLayout(sidebarPanel(id="sidebar",
                                                  basic.settings,
                                                  uiOutput("clonefreqs"),
                                                  actionButton("plus_button", "+")),
                                     mainPanel(plotOutput("path_plot"))
                       )
              ),
              tabPanel(title="Retrospective",
                       sidebarLayout(sidebarPanel(id="sidebar2",
                                                  basic.settings.2,
                                                  numericInput("cells_sequenced", "Number of cells sequenced", 1000, step=100, min=0),
                                                  uiOutput("clonefreqs2"),
                                                  actionButton("plus_button2", "+")),
                                     mainPanel(plotOutput("path_plot2"))
                       )
              ),
              tabPanel(title="FAQ", includeMarkdown("www/two_panel_faq.txt"))
  ),
  theme="theme.css"
)

server <- function(input, output, session)
{
  # Keep the basic settings in sync with each other across panels
  # When second cutoff is updated, if we're on the second tab, update the first cutoff
  observeEvent(input$cutoff2,
               {
                 if (input$tab == "Retrospective")
                 {
                   updateNumericInput(session, "cutoff", value=input$cutoff2)
                 }
               })
  # When the first cutoff is updated, if we're on the first tab, update the second tab
  observeEvent(input$cutoff,
               {
                 if (input$tab == "Prospective")
                 {
                   updateNumericInput(session, "cutoff2", value=input$cutoff)
                 }
               })
  # When second power is updated, if we're on the second tab, update the first cutoff
  observeEvent(input$power2,
               {
                 if (input$tab == "Retrospective")
                 {
                   updateNumericInput(session, "power", value=input$power2)
                 }
               })
  # When the first power is updated, if we're on the first tab, update the second tab
  observeEvent(input$power,
               {
                 if (input$tab == "Prospective")
                 {
                   updateNumericInput(session, "power2", value=input$power)
                 }
               })
  
  current.clones <- reactiveVal(1)
  current.clones.2 <- reactiveVal(1)
  
  # Keep track of which clone to remove
  clone.to.remove <- reactiveVal(NA)
  clone.to.remove.2 <- reactiveVal(NA)
  
  # Make an observer on the first minus button
  make.minus(1, input, session, clone.to.remove, current.clones, "minus", "clonefreq")
  
  # Same thing but for the second panel
  make.minus(1, input, session, clone.to.remove.2, current.clones.2, "retrominus", "retroclonefreq")
  
  # Respond to the plus button
  observeEvent(input$plus_button, {
    old.value <- current.clones()
    new.value <- old.value + 1
    current.clones(new.value)
    
    make.minus(new.value, input, session, clone.to.remove, current.clones, "minus", "clonefreq")
  })
  
  # Respond to the plus button on the second page
  observeEvent(input$plus_button2, {
    old.value <- current.clones.2()
    new.value <- old.value + 1
    current.clones.2(new.value)
    
    make.minus(new.value, input, session, clone.to.remove.2, current.clones.2, "retrominus", "retroclonefreq")
  })
  
  # Retrieve clone frequencies for first panel
  numeric.clonefreqs <- reactive({
    these.clonefreqs <- sapply(1:current.clones(), function(i) input[[sprintf("clonefreq%d",i)]])
    these.clonecounts <- sapply(1:current.clones(), function(i) input[[sprintf("clonecount%d",i)]])
    if (any(sapply(these.clonefreqs, is.null))) {
      NULL
    } else {
      rep.int(these.clonefreqs, these.clonecounts)
    }
  })
  # Same thing for the second panel. Don't need clone counts for this though, so it's a little different
  numeric.clonefreqs.2 <- reactive({
    these.clonefreqs <- sapply(1:current.clones.2(), function(i) input[[sprintf("retroclonefreq%d",i)]])
    if (any(sapply(these.clonefreqs, is.null))) {
      NULL
    } else {
      these.clonefreqs
    }
  })
  
  # Calculate the path. This was a bit much to copy and paste, so I made it a function.
  numeric.path <- numeric.path.maker(input, numeric.clonefreqs, current.clones)
  numeric.path.2 <- numeric.path.maker(input, numeric.clonefreqs.2, current.clones.2, "cells_sequenced")
  
  
  # Recommended numbers of cells to sequence for each panel
  recnum <- reactive({
    recnum <- which(numeric.path() >= input$power)[1] - 1
  })
  recnum.2 <- reactive({
    recnum.2 <- which(numeric.path.2() >= input$power2)[1] - 1
  })
  
  output$path_plot <- renderPlot(path.plot.maker(input, numeric.path, recnum, "power", "cutoff"))
  output$path_plot2 <- renderPlot(path.plot.maker(input, numeric.path.2, recnum.2, "power2", "cutoff2") + geom_vline(xintercept=input$cells_sequenced, colour="tomato2", size=rel(1.3)))
  
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
        if (i==1) {NULL} else {verticalLayout(actionButton(sprintf("minus%d", i), "-"))}
      )
    )
    
    # Since the clone has already been removed, set clone.to.remove to NA
    clone.to.remove(NA)
    
    # Put these boxes into a vertical layout and return
    do.call(verticalLayout, cancer.boxes)
  })
  
  
  clonefreq.boxes.2 <- eventReactive(current.clones.2(), {
    # What is the "current" clone frequency--that is, the frequency of the bottom clone on the list?
    clone.freqs <- lapply(1:current.clones.2(), function(i) input[[sprintf("retroclonefreq%d",i)]])
    if (all(sapply(clone.freqs, is.null))) {
      current <- .01
    } else {
      current <- clone.freqs[[max(which(!sapply(clone.freqs, is.null)))]]
    }
    
    # Make a mapping from new numbering to old numbering
    this.clone.to.remove <- clone.to.remove.2()
    if (is.na(this.clone.to.remove)) {
      old <- 1:current.clones.2()
    } else {
      before.removal <- 1:(current.clones.2()+1)
      old <- before.removal[-this.clone.to.remove]
    }
    # Make the boxes that allow user input of clone frequencies
    if (current.clones.2() == 0) {
      clone.indices <- integer(0)
    } else {
      clone.indices <- 1:current.clones.2()
    }
    cancer.boxes <- lapply(clone.indices, function(i)
      splitLayout(
        numericInput(sprintf("retroclonefreq%d", i), sprintf("Observed type %d", i),
                     ifelse(is.null(input[[sprintf("retroclonefreq%d",old[i])]]), current, input[[sprintf("retroclonefreq%d",old[i])]]), min=0.01, max=1, step=.01),
        actionButton(sprintf("retrominus%d", i), "-"),
        cellWidths=c("80%", "20%")
      )
    )
    
    # Since the clone has already been removed, set clone.to.remove to NA
    clone.to.remove.2(NA)
    
    # Put these boxes into a vertical layout and return
    do.call(verticalLayout, cancer.boxes)
  })
  
  
  output$clonefreqs <- renderUI({
    clonefreq.boxes()
  })
  output$clonefreqs2 <- renderUI({
    clonefreq.boxes.2()
  })
}

shinyApp(ui=ui, server=server)
