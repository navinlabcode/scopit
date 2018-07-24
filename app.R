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
  tabsetPanel(tabPanel(title="Calculation",
  sidebarLayout(sidebarPanel(id="sidebar",
                             numericInput("cutoff", "Number of cells of each type which must be sequenced", 3),
                             numericInput("power", "Required probability of sequencing this many cells from each type", .99,
                                          min=0, max=1, step=.01),
                             #numericInput("nclones", "Number of subpopulations of interest", 2, min=1),
                             uiOutput("clonefreqs"),
                             actionButton("plus_button", "+")),
#                             numericInput("normalfreq", "Remaining cells", NULL)),
                mainPanel(plotOutput("path_plot"))
  )
  ),
  tabPanel(title="Information",
           h4("What is this tool used for?"),
           p("This tool is for estimating how many cells must be sequenced in a single-cell sequencing experiment. It is intended for rough 'napkin calculations', since it is based on information about subpopulation frequencies that you probably don't really have."),
           h4("Who made it?"),
           HTML("<p>Alexander Davis, a Ph.D. candidate in the laboratory of Dr. Nicholas Navin, at the University of Texas MD Anderson Cancer Center Department of Genetics. If you have questions, bug reports, or feature requests, contact me at: <img src='image.png'>.</p>"),
           h4("Is it for single-cell DNA sequencing or single-cell RNA sequencing?"),
           p("Either."),
           h4("What constitutes a 'subpopulation'?"),
           p("In cancer genomics, a 'subpopulation' would be a subclone of a tumor. For experiments on normal tissue, a 'subpopulation' is typically a cell type that the researcher wishes to discover, and the 'remaining cells' would be known or otherwise uninteresting cell types."),
           h4("If I set the number of subpopulations to n, why do I see n+1 boxes below it?"),
           p(HTML("You set the number of subpopulations <i>of interest</i> to n. The last box is the frequency of cells which are not of interest. In calculating the number of cells required, it is assumed that no cells are required from these remaining cells.")),
           h4("How many cells do I want to sequence from each subpopulation?"),
           p("It should be the size of the smallest cluster you can detect in your clustering analysis. In single-cell DNA sequencing for cancer genomics, it may be just a few, whereas for RNA sequencing, it may be much more."),
           h4("What should I set the probability to?"),
           p("Since this is the probability that your very expensive experiment meets its stated goals, I would recommend something very high, like 99%. However, in the analogous setting of power calculation, 80% is also a common choice."),
           h4("How do I know how many subpopulations there are, and what their frequencies are?"),
           p("You probably don't, or you wouldn't have to do the experiment. However, with this tool you can experiment with different plausible frequencies to understand what they imply about the required number of cells."),
           h4("Why can't I enter frequencies below 0.001?"),
           p("Because the tool is currently not fast enough to handle them. Increasing speed is a priority for future versions. I expect that I will have to start using an approximate method, but I think with the right approximation it is possible to guarantee at least several digits of precision."),
           h4("What is the tool actually doing?"),
           p("It is calculating a probability in a multinomial distribution. Specifically, it is calculating"),
           withMathJax(p("$$P(X_1 \\ge c, \\dots, X_k \\ge c, \\Omega \\ge 0)$$")),
           withMathJax(p("where \\(X_i\\) is the number of cells observed from subpopulation \\(i\\), \\(\\Omega\\) is the number of cells observed outside of these subpopulations, and \\(c\\) is the number of cells that must be sequenced from each subpopulation. \\(X_1, \\dots,X_k,\\Omega\\) jointly have a multinomial distribution with \\(n\\) equal to the sample size, and probabilities equal to the specified subpopulation frequencies.")),
           p(HTML("The calculation of this multinomial probability follows a strategy suggested by <a href='http://dx.doi.org/10.1214/aos/1176345593'>Bruce Levin (1981) <i>Annals of Statistics</i></a>. First, the probability of the same event is calculated under the assumption that each number is independently Poisson distributed. Then, this probability is updated to a probability conditional on the total number of cells, using Bayes' theorem. Since Poisson random variables conditional on their sum have a multinomial distribution, this gives the desired probability. The bottleneck step of the calculation is the likelihood in Bayes' theorem, which is a probability of a sum of truncated Poisson distributions. Levin approximates this probability, but suggests that an exact calculation could be accomplished using convolutions, which is the strategy used here. Convolutions are calculated using the fast Fourier transform, using the R function 'convolve'.")),
           h4("Is the resulting number an approximation?"),
           p("No, it is exact. Any number calculated by a computer has small round-off errors, but no approximation has been made in the math."),
           h4("How can I save my query and the results?"),
           p("You currently can't, except by right-clicking the plot and saving it the way you would save any other image. This feature is a priority for future versions."),
           h4("Where can I find examples of these kinds of calculations?"),
           p(HTML("This kind of analysis can be found in Supplementary Figure 9 of <a href='http://dx.doi.org/10.1038/ng.3641'>Gao et al. (2016) <i>Nature Genetics</i></a> and Supplementary Figure 7 of <a href='https://doi.org/10.1016/j.cell.2017.12.007'>Casasent et al. (2018) <i>Cell</i></a>.")),
           h4("Why do the horizontal and vertical lines meet a little bit below the curve?"),
           p("The probability you input may not be exactly achievable, since the number of cells must be a whole number. If this is the case, then this tool gives you the smallest sample size at which the probability exceeds your input. This means the curve (the probability) is above the horizontal line (your input)."),
           h4("Where can I find the code?"),
           p(HTML("The git repository for this app is hosted at <a href='https://github.com/alexdavisscs/shiny-multinomial'>https://github.com/alexdavisscs/shiny-multinomial</a>.")),
           h4("Acknowledgements"),
           p("The idea of a Shiny app for multinomial calculations is due to Dr. Ruli Gao, who also contributed her own calculations to check the accuracy of the algorithm. Aislyn Schalck gave extensive comments and advice. Charissa Kim tested early versions and provided encouragement. Dr. Nicholas Navin has guided all of my work in single cell sequencing. I (Alexander Davis) am supported by the National Library of Medicine Training Program in Biomedical Informatics (4T15LM007093-25), as well as the American Legion Auxiliary. Dr. Navin is an Andrew Sabin Family Fellow. This research was also supported by grants to Dr. Navin from the National Cancer Institute (1RO1CA169244-01) and American Cancer Society (129098-RSG-16-092-01-TBG)."))),
                theme="theme.css"
)

server <- function(input, output, session)
{
  current.clones <- reactiveVal(1)
  
  # Keep track of which clone to remove
  clone.to.remove <- reactiveVal(NA)
  
  # Make an observer on the first minus button
  make.minus(1, input, session, clone.to.remove, current.clones)

  # Old code to respond to the number of clones
#  observeEvent(input$nclones, {
#    if (!is.na(input$nclones))
#    {
#      current.clones(input$nclones)
#    }
#  })
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
    #if (any(numeric.clonefreqs() < .001)) stop("Minimum subpopulation frequency is .001 until a faster algorithm is implemented")
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
      current <- .1
    } else {
      current <- clone.freqs[[max(which(!sapply(clone.freqs, is.null)))]]
      #current <- clone.freqs[[current.clones()]]
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
      splitLayout(
        numericInput(sprintf("clonefreq%d", i), sprintf("Frequency %d", i),
                     ifelse(is.null(input[[sprintf("clonefreq%d",old[i])]]), current, input[[sprintf("clonefreq%d",old[i])]]), min=0.01, max=1, step=.01),
        numericInput(sprintf("clonecount%d", i), sprintf("# of types", i), 1),
        actionButton(sprintf("minus%d", i), "-")
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

  # Enforce the value of normalfreq when the user tries to change it
  observeEvent(input$normalfreq, {
    normal.cell.freq <- 1 - sum(numeric.clonefreqs())
    updateNumericInput(session, "normalfreq", value=normal.cell.freq)
  })

  # Keep normal cell frequency up to date when clones are changed
  observe({
    normal.cell.freq <- 1 - sum(numeric.clonefreqs())
    updateNumericInput(session, "normalfreq", value=normal.cell.freq)
  })
}

shinyApp(ui=ui, server=server)
