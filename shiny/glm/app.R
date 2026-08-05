# ---------------------------------------------------------------------------
# GLM simulator  --  SBOEM112, part "Linear statistical modeling"
# F. De Laender, UNamur
#
# Simulates data from a GLM using the notation of the course:
#
#     y_i ~ D                          (random component)
#     eta_i = sum_j beta_j x_ij        (linear predictor)
#     g(mu_i) = eta_i                  (canonical link function)
#
# SIMULATION ONLY. No model is fitted anywhere in this app. The curve drawn
# on the plot is the TRUE mean mu_i implied by the parameters you chose --
# it is not an estimate.
# ---------------------------------------------------------------------------

library(shiny)
library(ggplot2)
library(DT)

# ---------------------------------------------------------------------------
# 1. SIMULATION ENGINE
#
#    Everything that is specific to a distribution D lives in DISTRIBUTIONS,
#    in the same order as the tables in the slides (tbl-theta and the
#    canonical link table). To add a distribution, add one entry.
#
#    To extend to LMMs later: give simulate_glm() a Z matrix and a vector u
#    of random effects, and change one line -- eta <- X %*% beta + Z %*% u.
#    Nothing else in this file needs to know about it.
# ---------------------------------------------------------------------------

DISTRIBUTIONS <- list(
  normal = list(
    label      = "Normal",
    link_label = "identity",
    link_math  = "\\mu_i = \\eta_i",
    var_math   = "\\sigma^2",
    linkinv    = function(eta, pars) eta,
    variance   = function(mu, pars) rep(pars$sigma^2, length(mu)),
    rand       = function(mu, pars) rnorm(length(mu), mean = mu, sd = pars$sigma),
    uses_sigma = TRUE,
    uses_n     = FALSE
  ),
  poisson = list(
    label      = "Poisson",
    link_label = "log",
    link_math  = "\\log(\\mu_i) = \\eta_i",
    var_math   = "\\mu_i",
    linkinv    = function(eta, pars) exp(eta),
    variance   = function(mu, pars) mu,
    rand       = function(mu, pars) rpois(length(mu), lambda = mu),
    uses_sigma = FALSE,
    uses_n     = FALSE
  ),
  binom_prop = list(
    label      = "Binomial (proportion)",
    link_label = "logit",
    link_math  = "\\log\\left(\\frac{\\mu_i}{1-\\mu_i}\\right) = \\eta_i",
    var_math   = "\\mu_i(1-\\mu_i)/n_i",
    linkinv    = function(eta, pars) plogis(eta),
    variance   = function(mu, pars) mu * (1 - mu) / pars$n,
    rand       = function(mu, pars) rbinom(length(mu), size = pars$n, prob = mu) / pars$n,
    uses_sigma = FALSE,
    uses_n     = TRUE
  ),
  binom_count = list(
    label      = "Binomial (nr of successes)",
    link_label = "logit",
    link_math  = "\\log\\left(\\frac{\\pi_i}{1-\\pi_i}\\right) = \\eta_i",
    var_math   = "n_i\\pi_i(1-\\pi_i)",
    linkinv    = function(eta, pars) pars$n * plogis(eta),
    variance   = function(mu, pars) mu * (1 - mu / pars$n),
    rand       = function(mu, pars) rbinom(length(mu), size = pars$n, prob = mu / pars$n),
    uses_sigma = FALSE,
    uses_n     = TRUE
  )
)

DIST_CHOICES <- setNames(names(DISTRIBUTIONS), vapply(DISTRIBUTIONS, `[[`, "", "label"))

# Build the design: distinct x values x (optional) 2-level factor x replicates.
# Columns are named with the course notation: x_i1 = 1 (intercept), x_i2 = the
# continuous predictor, x_i3 = indicator for level B, x_i4 = the interaction.
build_design <- function(x_values, use_factor, n_rep) {
  groups <- if (use_factor) c("A", "B") else "A"
  d <- expand.grid(
    replicate = seq_len(n_rep),
    group     = groups,
    x         = x_values,
    KEEP.OUT.ATTRS   = FALSE,
    stringsAsFactors = FALSE
  )
  d <- d[order(d$x, d$group, d$replicate), , drop = FALSE]
  d$x_i1 <- 1
  d$x_i2 <- d$x
  d$x_i3 <- as.numeric(d$group == "B")
  d$x_i4 <- d$x_i2 * d$x_i3
  rownames(d) <- NULL
  d
}

simulate_glm <- function(design, beta, dist_key, pars) {
  spec <- DISTRIBUTIONS[[dist_key]]
  X   <- as.matrix(design[, c("x_i1", "x_i2", "x_i3", "x_i4"), drop = FALSE])
  eta <- as.numeric(X %*% beta)
  mu  <- spec$linkinv(eta, pars)
  v   <- spec$variance(mu, pars)
  y   <- spec$rand(mu, pars)
  data.frame(
    i     = seq_len(nrow(design)),
    x_i1  = design$x_i1,
    x_i2  = design$x_i2,
    x_i3  = design$x_i3,
    x_i4  = design$x_i4,
    group = design$group,
    eta_i = eta,
    mu_i  = mu,
    var_i = v,
    y_i   = y,
    stringsAsFactors = FALSE
  )
}

# The true mean curve, on a fine grid, for drawing on the plot.
mean_curve <- function(x_grid, beta, dist_key, pars, use_factor) {
  groups <- if (use_factor) c("A", "B") else "A"
  g <- expand.grid(x_i2 = x_grid, group = groups,
                   KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  ind <- as.numeric(g$group == "B")
  X   <- cbind(1, g$x_i2, ind, g$x_i2 * ind)
  g$mu_i <- DISTRIBUTIONS[[dist_key]]$linkinv(as.numeric(X %*% beta), pars)
  g
}

GROUP_COLS <- c(A = "#0072B2", B = "#D55E00")

# ---------------------------------------------------------------------------
# 2. USER INTERFACE
# ---------------------------------------------------------------------------

ui <- fluidPage(
  withMathJax(),
  titlePanel("Simulating data from a GLM"),
  tags$p(
    style = "color:#555; margin-top:-8px;",
    "SBOEM112 – Linear statistical modeling. ",
    tags$em("No model is fitted here: you choose the parameters, the app generates the data.")
  ),
  sidebarLayout(
    sidebarPanel(
      width = 4,

      tags$h4("1. Random component"),
      selectInput("dist", "Distribution \\(D\\)", choices = DIST_CHOICES),
      conditionalPanel(
        "input.dist == 'normal'",
        sliderInput("sigma", "\\(\\sigma\\) (residual sd)", min = 0.05, max = 5,
                    value = 1, step = 0.05)
      ),
      conditionalPanel(
        "input.dist == 'binom_prop' || input.dist == 'binom_count'",
        sliderInput("n_trials", "\\(n_i\\) (nr of trials)", min = 1, max = 50,
                    value = 10, step = 1)
      ),

      tags$hr(),
      tags$h4("2. Predictors"),
      sliderInput("n_x", "Nr of distinct values of \\(x_{i2}\\)",
                  min = 1, max = 10, value = 6, step = 1),
      sliderInput("x_range", "Range of \\(x_{i2}\\)",
                  min = -5, max = 10, value = c(0, 5), step = 0.5),
      sliderInput("n_rep", "Replicates per predictor value",
                  min = 1, max = 50, value = 10, step = 1),
      checkboxInput("use_factor", "Add a 2-level factor (levels A and B)", value = FALSE),
      conditionalPanel(
        "input.use_factor == true",
        checkboxInput("use_interaction", "Let the factor change the slope (interaction)",
                      value = FALSE)
      ),

      tags$hr(),
      tags$h4("3. Free parameters"),
      sliderInput("beta1", "\\(\\beta_1\\) (intercept)", min = -5, max = 5,
                  value = 0.5, step = 0.1),
      sliderInput("beta2", "\\(\\beta_2\\) (effect of \\(x_{i2}\\))", min = -5, max = 5,
                  value = 0.4, step = 0.1),
      conditionalPanel(
        "input.use_factor == true",
        sliderInput("beta3", "\\(\\beta_3\\) (effect of level B)", min = -5, max = 5,
                    value = 1, step = 0.1)
      ),
      conditionalPanel(
        "input.use_factor == true && input.use_interaction == true",
        sliderInput("beta4", "\\(\\beta_4\\) (interaction)", min = -5, max = 5,
                    value = 0, step = 0.1)
      ),

      tags$hr(),
      actionButton("resample", "Draw a new sample", class = "btn-primary"),
      tags$p(style = "color:#777; font-size:90%; margin-top:6px;",
             "Same parameters, new random draw.")
    ),

    mainPanel(
      width = 8,
      uiOutput("model_box"),
      tabsetPanel(
        tabPanel(
          "Plot",
          br(),
          plotOutput("plot", height = "460px"),
          tags$p(style = "color:#777; font-size:90%;",
                 "Points are the simulated \\(y_i\\) (jittered horizontally so replicates ",
                 "are visible). The line is the true mean \\(\\mu_i\\).")
        ),
        tabPanel(
          "Simulated data",
          br(),
          DT::dataTableOutput("table")
        )
      )
    )
  )
)

# ---------------------------------------------------------------------------
# 3. SERVER
# ---------------------------------------------------------------------------

server <- function(input, output, session) {

  seed <- reactiveVal(1)
  observeEvent(input$resample, seed(seed() + 1))

  spec <- reactive(DISTRIBUTIONS[[input$dist]])

  pars <- reactive(list(
    sigma = if (is.null(input$sigma)) 1 else input$sigma,
    n     = if (is.null(input$n_trials)) 10 else input$n_trials
  ))

  beta <- reactive({
    use_f <- isTRUE(input$use_factor)
    use_i <- use_f && isTRUE(input$use_interaction)
    b3 <- if (use_f && !is.null(input$beta3)) input$beta3 else 0
    b4 <- if (use_i && !is.null(input$beta4)) input$beta4 else 0
    c(input$beta1, input$beta2, b3, b4)
  })

  x_values <- reactive({
    r <- input$x_range
    if (input$n_x == 1) r[1] else seq(r[1], r[2], length.out = input$n_x)
  })

  sim <- reactive({
    d <- build_design(x_values(), isTRUE(input$use_factor), input$n_rep)
    set.seed(seed())
    out <- simulate_glm(d, beta(), input$dist, pars())
    validate(need(
      all(is.finite(out$mu_i)) && max(out$mu_i) < 1e6,
      paste("The mean mu_i has blown up. Remember that for this distribution the",
            "betas act on the link scale, so eta_i is exponentiated. Lower beta1",
            "or beta2, or narrow the range of x.")
    ))
    out
  })

  curve_df <- reactive({
    xv <- x_values()
    grid <- if (length(xv) == 1) xv else seq(min(xv), max(xv), length.out = 200)
    mean_curve(grid, beta(), input$dist, pars(), isTRUE(input$use_factor))
  })

  # -- the model, written in course notation, updating with the choices ------
  output$model_box <- renderUI({
    use_f <- isTRUE(input$use_factor)
    use_i <- use_f && isTRUE(input$use_interaction)
    eta_terms <- "\\beta_1 x_{i1} + \\beta_2 x_{i2}"
    if (use_f) eta_terms <- paste0(eta_terms, " + \\beta_3 x_{i3}")
    if (use_i) eta_terms <- paste0(eta_terms, " + \\beta_4 x_{i4}")

    key <- paste0(
      "\\(x_{i1}=1\\), \\(x_{i2}\\) = the continuous predictor",
      if (use_f) ", \\(x_{i3}=1\\) for level B and \\(0\\) for level A" else "",
      if (use_i) ", \\(x_{i4}=x_{i2}x_{i3}\\)" else ""
    )

    withMathJax(div(
      style = "background:#f7f7f7; border-left:4px solid #0072B2; padding:10px 14px; margin-bottom:14px;",
      div(sprintf("$$y_i \\sim \\text{%s} \\\\ \\eta_i = %s \\\\ %s$$",
                  spec()$label, eta_terms, spec()$link_math)),
      div(style = "font-size:90%; color:#555;",
          HTML(paste0("with ", key, ".")),
          br(),
          sprintf("Canonical link: %s. Variance: \\(\\text{var}(y_i) = %s\\).",
                  spec()$link_label, spec()$var_math))
    ))
  })

  output$plot <- renderPlot({
    d  <- sim()
    cd <- curve_df()
    use_f <- isTRUE(input$use_factor)

    xr <- diff(range(d$x_i2))
    jw <- if (xr > 0 && input$n_x > 1) (xr / (input$n_x - 1)) * 0.12 else 0.05

    p <- ggplot() +
      theme_bw(base_size = 14) +
      labs(x = expression(x[i2]), y = expression(y[i]), colour = "group")

    if (use_f) {
      p <- p +
        geom_point(data = d, aes(x = x_i2, y = y_i, colour = group),
                   position = position_jitter(width = jw, height = 0, seed = 1),
                   alpha = 0.65, size = 2) +
        scale_colour_manual(values = GROUP_COLS)
    } else {
      p <- p +
        geom_point(data = d, aes(x = x_i2, y = y_i),
                   position = position_jitter(width = jw, height = 0, seed = 1),
                   alpha = 0.65, size = 2, colour = GROUP_COLS[["A"]])
    }

    if (nrow(cd) > (if (use_f) 2 else 1)) {
      p <- if (use_f) {
        p + geom_line(data = cd, aes(x = x_i2, y = mu_i, colour = group), linewidth = 1.1)
      } else {
        p + geom_line(data = cd, aes(x = x_i2, y = mu_i), linewidth = 1.1,
                      colour = GROUP_COLS[["A"]])
      }
    } else {
      p <- if (use_f) {
        p + geom_point(data = cd, aes(x = x_i2, y = mu_i, colour = group),
                       size = 6, shape = 95)
      } else {
        p + geom_point(data = cd, aes(x = x_i2, y = mu_i), size = 6, shape = 95,
                       colour = GROUP_COLS[["A"]])
      }
    }
    p
  })

  output$table <- DT::renderDataTable({
    d <- sim()
    keep <- c("i", "x_i1", "x_i2")
    if (isTRUE(input$use_factor)) keep <- c(keep, "group", "x_i3", "x_i4")
    keep <- c(keep, "eta_i", "mu_i", "var_i", "y_i")
    d <- d[, keep, drop = FALSE]
    names(d)[names(d) == "var_i"] <- "var(y_i)"

    DT::datatable(
      d,
      rownames = FALSE,
      options  = list(
        paging = FALSE, searching = FALSE, info = FALSE,
        scrollY = "420px", scrollX = TRUE, scrollCollapse = TRUE,
        dom = "t"
      )
    ) |>
      DT::formatRound(columns = c("x_i2", "eta_i", "mu_i", "var(y_i)", "y_i"), digits = 3)
  })
}

shinyApp(ui, server)
