library(shiny)

shinyUI(fluidPage(
  titlePanel("Ivy: \"Incremental Validity\" Error Rate Calculator"),
  p("Westfall, J., & Yarkoni, T. (2016). Statistically controlling for confounding constructs is harder than you think. ",em("PLoS ONE, 11"),"(3), e0152719."),
  fluidRow(
    column(1, div(a("Article", href="http://dx.doi.org/10.1371/journal.pone.0152719"))),
    column(3, div(a("Supplemental Appendix S1: Derivation of statistical properties of incremental validity", href="http://jakewestfall.org/publications/ivy_appendix.pdf"))),
    column(2, div(a("Code for this app (using package 'shiny' in R)", href="https://github.com/jake-westfall/ivy_app"))),
    column(2, div(a("Back to JakeWestfall.org", href="http://jakewestfall.org")))
  ),
  helpText("Note: when sharing the link to this app, please use the stable redirecting page at",
           a("jakewestfall.org/ivy/,", href="http://jakewestfall.org/ivy/"),
           "as the app's current URL is possibly subject to change."),

  sidebarLayout(
    
    # INPUT
    sidebarPanel(
      numericInput("rho1", label="rho1 (Y-T1 correlation)",
                   value=.5, min=-1, max=1, step=.01),
      radioButtons("rho1_type", label=NULL, inline=TRUE,
                   choices=list("Simple"=1,"Partial"=2)),
      numericInput("rho2", label="rho2 (Y-T2 correlation)",
                   value=.5, min=-1, max=1, step=.01),
      radioButtons("rho2_type", label=NULL, inline=TRUE,
                   choices=list("Simple"=1,"Partial"=2)),
      numericInput("delta", label="delta (T1-T2 correlation)",
                   value=.5, min=-1, max=1, step=.01),
      radioButtons("delta_type", label=NULL, inline=TRUE,
                   choices=list("Simple"=1,"Partial"=2)),
      numericInput("alpha1", label="alpha1 (Reliabiity of X1)",
                   value=.5, min=0, max=1, step=.01),
      numericInput("alpha2", label="alpha2 (Reliability of X2)",
                   value=.5, min=0, max=1, step=.01),
      numericInput("n", label="n (Sample size)", value=100, min=10, step=10)
    ),
    
    # OUTPUT
    mainPanel(
      htmlOutput("error"),
      
      h3("Probabilities of rejecting regression coefficients"),
      strong("Probability of rejecting both Beta1=0 and Beta2=0:"),
      textOutput("pBoth"),
      strong("Probability of rejecting only Beta1=0:"),
      textOutput("pB1"),
      strong("Probability of rejecting only Beta2=0:"),
      textOutput("pB2"),
      strong("Probability of rejecting neither:"),
      textOutput("pNone"),
      
      fluidRow(
        column(5, h3("Latent variable parameters")),
        column(5, h3("Observed variable parameters"))
      ),
      fluidRow(
        column(5, strong("rho1_simple:")),
        column(5, strong("Y-X1 simple correlation:"))
      ),
      fluidRow(
        column(5, textOutput("rho1_simple")),
        column(5, textOutput("scorr1"))
      ),
      fluidRow(
        column(5, strong("rho1_partial:")),
        column(5, strong("Y-X1 partial correlation:"))
      ),
      fluidRow(
        column(5, textOutput("rho1_partial")),
        column(5, textOutput("pcorr1"))
      ),
      fluidRow(
        column(5, strong("rho2_simple:")),
        column(5, strong("Y-X2 simple correlation:"))
      ),
      fluidRow(
        column(5, textOutput("rho2_simple")),
        column(5, textOutput("scorr2"))
      ),
      fluidRow(
        column(5, strong("rho2_partial:")),
        column(5, strong("Y-X2 partial correlation:"))
      ),
      fluidRow(
        column(5, textOutput("rho2_partial")),
        column(5, textOutput("pcorr2"))
      ),
      fluidRow(
        column(5, strong("delta_simple:")),
        column(5, strong("X1-X2 simple correlation:"))
      ),
      fluidRow(
        column(5, textOutput("delta_simple")),
        column(5, textOutput("r_x1x2"))
      ),
      fluidRow(
        column(5, strong("delta_partial:")),
        column(5, strong("X1-X2 partial correlation:"))
      ),
      fluidRow(
        column(5, textOutput("delta_partial")),
        column(5, textOutput("pcorr3"))
      )
    ) # end mainPanel
  ) # end sidebarLayout
)) # end shinyUI and fluidPage