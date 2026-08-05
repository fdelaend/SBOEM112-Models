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
# SIMULATION ONLY. No model is fitted anywhere in this app. The distributions
# drawn on the plot are the TRUE distributions D implied by the parameters you
# chose -- they are not estimates.
# ---------------------------------------------------------------------------

library(shiny)
library(ggplot2)
library(DT)

BLOWUP_MSG <- paste(
  "The mean mu_i has blown up. Remember that for this distribution the betas",
  "act on the link scale, so eta_i gets exponentiated. Lower beta1 or beta2,",
  "or narrow the range of x."
)

# ---------------------------------------------------------------------------
# 1. SIMULATION ENGINE
#
#    Everything that is specific to a distribution D lives in DISTRIBUTIONS,
#    in the same order as the tables in the slides (tbl-theta and the
#    canonical link table). To add a distribution, add one entry:
#
#      linkinv   inverse of the canonical link: eta -> mu
#      variance  var(y_i) as a function of mu
#      rand      random draw from D
#      dens      the distribution itself: returns data.frame(y, d), where d is
#                a density (continuous D) or a probability (discrete D)
#      discrete  TRUE -> drawn as spikes, FALSE -> drawn as a smooth curve
#
#    To extend to LMMs later, change one line in simulate_glm():
#    eta <- X %*% beta + Z %*% u.
# ---------------------------------------------------------------------------

DISTRIBUTIONS <- list(
  normal = list(
    label      = "Normal",
    link_label = "identity",
    link_math  = "\\mu_i = \\eta_i",
    var_math   = "\\sigma^2",
    discrete   = FALSE,
    linkinv    = function(eta, pars) eta,
    variance   = function(mu, pars) rep(pars$sigma^2, length(mu)),
    rand       = function(mu, pars) rnorm(length(mu), mean = mu, sd = pars$sigma),
    dens       = function(mu, pars) {
      y <- seq(mu - 4 * pars$sigma, mu + 4 * pars$sigma, length.out = 200)
      data.frame(y = y, d = dnorm(y, mean = mu, sd = pars$sigma))
    },
    uses_sigma = TRUE,
    uses_n     = FALSE
  ),
  poisson = list(
    label      = "Poisson",
    link_label = "log",
    link_math  = "\\log(\\mu_i) = \\eta_i",
    var_math   = "\\mu_i",
    discrete   = TRUE,
    linkinv    = function(eta, pars) exp(eta),
    variance   = function(mu, pars) mu,
    rand       = function(mu, pars) rpois(length(mu), lambda = mu),
    dens       = function(mu, pars) {
      top <- max(1, min(qpois(0.9999, mu), 500))
      y   <- 0:top
      data.frame(y = y, d = dpois(y, lambda = mu))
    },
    uses_sigma = FALSE,
    uses_n     = FALSE
  ),
  binom_prop = list(
    label      = "Binomial (proportion)",
    link_label = "logit",
    link_math  = "\\log\\left(\\frac{\\mu_i}{1-\\mu_i}\\right) = \\eta_i",
    var_math   = "\\mu_i(1-\\mu_i)/n_i",
    discrete   = TRUE,
    linkinv    = function(eta, pars) plogis(eta),
    variance   = function(mu, pars) mu * (1 - mu) / pars$n,
    rand       = function(mu, pars) rbinom(length(mu), size = pars$n, prob = mu) / pars$n,
    dens       = function(mu, pars) {
      k <- 0:pars$n
      data.frame(y = k / pars$n, d = dbinom(k, size = pars$n, prob = mu))
    },
    uses_sigma = FALSE,
    uses_n     = TRUE
  ),
  binom_count = list(
    label      = "Binomial (nr of successes)",
    link_label = "logit",
    link_math  = "\\log\\left(\\frac{\\pi_i}{1-\\pi_i}\\right) = \\eta_i",
    var_math   = "n_i\\pi_i(1-\\pi_i)",
    discrete   = TRUE,
    linkinv    = function(eta, pars) pars$n * plogis(eta),
    variance   = function(mu, pars) mu * (1 - mu / pars$n),
    rand       = function(mu, pars) rbinom(length(mu), size = pars$n, prob = mu / pars$n),
    dens       = function(mu, pars) {
      k <- 0:pars$n
      data.frame(y = k, d = dbinom(k, size = pars$n, prob = mu / pars$n))
    },
    uses_sigma = FALSE,
    uses_n     = TRUE
  )
)

DIST_CHOICES <- setNames(names(DISTRIBUTIONS), vapply(DISTRIBUTIONS, `[[`, "", "label"))
GROUP_COLS   <- c(A = "#0072B2", B = "#D55E00")

# Build the design: distinct x values x (optional) 2-level factor x replicates.
# Columns use the course notation: x_i1 = 1 (intercept), x_i2 = the continuous
# predictor, x_i3 = indicator for level B, x_i4 = the interaction.
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
  data.frame(
    i     = seq_len(nrow(design)),
    x_i1  = design$x_i1,
    x_i2  = design$x_i2,
    x_i3  = design$x_i3,
    x_i4  = design$x_i4,
    group = design$group,
    eta_i = eta,
    mu_i  = mu,
    var_i = spec$variance(mu, pars),
    y_i   = spec$rand(mu, pars),
    stringsAsFactors = FALSE
  )
}

# Pick k roughly evenly spaced values out of xv.
pick_x <- function(xv, k) {
  k <- min(k, length(xv))
  if (k < 1) return(numeric(0))
  xv[unique(round(seq(1, length(xv), length.out = k)))]
}

# The true distribution of y_i at selected x values, ready to be drawn
# sideways: each shape runs from its baseline at x0 out to x0 + width.
distribution_shapes <- function(x_at, beta, dist_key, pars, use_factor, spacing) {
  spec   <- DISTRIBUTIONS[[dist_key]]
  groups <- if (use_factor) c("A", "B") else "A"
  out    <- list()
  for (g in groups) {
    ind <- as.numeric(g == "B")
    for (x0 in x_at) {
      eta <- sum(beta * c(1, x0, ind, x0 * ind))
      mu  <- spec$linkinv(eta, pars)
      if (!is.finite(mu) || mu > 1e6) return(NULL)
      dd <- spec$dens(mu, pars)
      out[[length(out) + 1L]] <- data.frame(
        x0 = x0, group = g, mu = mu, y = dd$y, d = dd$d,
        stringsAsFactors = FALSE
      )
    }
  }
  res <- do.call(rbind, out)
  # One common scale for all shapes, so their relative heights stay
  # comparable -- this is what makes the change in the distribution visible.
  maxd <- max(res$d, na.rm = TRUE)
  if (!is.finite(maxd) || maxd <= 0) maxd <- 1
  res$xend <- res$x0 + (res$d / maxd) * spacing * 0.85
  res
}

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
          wellPanel(
            style = "padding:10px 14px;",
            radioButtons("show", "Show",
                         choices  = c("Simulated data"     = "data",
                                      "True distributions" = "dens",
                                      "Both"               = "both"),
                         selected = "both", inline = TRUE),
            conditionalPanel(
              "input.show != 'data'",
              sliderInput("n_dens", "At how many values of \\(x_{i2}\\)?",
                          min = 1, max = 10, value = 4, step = 1, width = "320px")
            )
          ),
          plotOutput("plot", height = "460px"),
          uiOutput("plot_note")
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

  # Spacing between neighbouring x values: sets how wide a distribution may be.
  spacing <- reactive({
    xv <- x_values()
    if (length(xv) > 1) min(diff(sort(xv))) else 1
  })

  sim <- reactive({
    d <- build_design(x_values(), isTRUE(input$use_factor), input$n_rep)
    set.seed(seed())
    out <- simulate_glm(d, beta(), input$dist, pars())
    validate(need(all(is.finite(out$mu_i)) && max(out$mu_i) < 1e6, BLOWUP_MSG))
    out
  })

  shapes <- reactive({
    x_at <- pick_x(x_values(), input$n_dens)
    res  <- distribution_shapes(x_at, beta(), input$dist, pars(),
                                isTRUE(input$use_factor), spacing())
    validate(need(!is.null(res), BLOWUP_MSG))
    res
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

  output$plot_note <- renderUI({
    bits <- character(0)
    if (input$show %in% c("data", "both"))
      bits <- c(bits, paste("Points are the simulated \\(y_i\\), jittered horizontally",
                            "so replicates are visible."))
    if (input$show %in% c("dens", "both"))
      bits <- c(bits, paste(
        "Each sideways shape is the true distribution of \\(y_i\\) at that value of",
        "\\(x_{i2}\\), drawn from its baseline to the right. All shapes share one",
        "scale, so their widths are comparable. The dot marks the mean \\(\\mu_i\\)."))
    withMathJax(tags$p(style = "color:#777; font-size:90%;", paste(bits, collapse = " ")))
  })

  output$plot <- renderPlot({
    use_f    <- isTRUE(input$use_factor)
    show_pts <- input$show %in% c("data", "both")
    show_den <- input$show %in% c("dens", "both")

    p <- ggplot() +
      theme_bw(base_size = 14) +
      labs(x = expression(x[i2]), y = expression(y[i]), colour = "group")

    # -- the true distributions, drawn sideways -----------------------------
    if (show_den) {
      sh <- shapes()
      base <- do.call(rbind, lapply(
        split(sh, list(sh$x0, sh$group), drop = TRUE),
        function(z) data.frame(x0 = z$x0[1], group = z$group[1], mu = z$mu[1],
                               ymin = min(z$y), ymax = max(z$y),
                               stringsAsFactors = FALSE)
      ))

      # baseline of each distribution
      p <- if (use_f) {
        p + geom_segment(data = base,
                         aes(x = x0, xend = x0, y = ymin, yend = ymax, colour = group),
                         linewidth = 0.3, alpha = 0.5)
      } else {
        p + geom_segment(data = base,
                         aes(x = x0, xend = x0, y = ymin, yend = ymax),
                         linewidth = 0.3, alpha = 0.5, colour = GROUP_COLS[["A"]])
      }

      # the distribution itself: spikes if D is discrete, a curve if continuous
      if (isTRUE(spec()$discrete)) {
        p <- if (use_f) {
          p + geom_segment(data = sh,
                           aes(x = x0, xend = xend, y = y, yend = y, colour = group),
                           linewidth = 0.9, alpha = 0.85)
        } else {
          p + geom_segment(data = sh,
                           aes(x = x0, xend = xend, y = y, yend = y),
                           linewidth = 0.9, alpha = 0.85, colour = GROUP_COLS[["A"]])
        }
      } else {
        p <- if (use_f) {
          p + geom_path(data = sh,
                        aes(x = xend, y = y, colour = group,
                            group = interaction(x0, group)),
                        linewidth = 0.9, alpha = 0.9)
        } else {
          p + geom_path(data = sh,
                        aes(x = xend, y = y, group = x0),
                        linewidth = 0.9, alpha = 0.9, colour = GROUP_COLS[["A"]])
        }
      }

      # the mean
      p <- if (use_f) {
        p + geom_point(data = base, aes(x = x0, y = mu, colour = group), size = 2.4)
      } else {
        p + geom_point(data = base, aes(x = x0, y = mu), size = 2.4,
                       colour = GROUP_COLS[["A"]])
      }
    }

    # -- the simulated data -------------------------------------------------
    if (show_pts) {
      d  <- sim()
      jw <- if (input$n_x > 1) spacing() * 0.12 else 0.05
      p <- if (use_f) {
        p + geom_point(data = d, aes(x = x_i2, y = y_i, colour = group),
                       position = position_jitter(width = jw, height = 0, seed = 1),
                       alpha = 0.55, size = 2)
      } else {
        p + geom_point(data = d, aes(x = x_i2, y = y_i),
                       position = position_jitter(width = jw, height = 0, seed = 1),
                       alpha = 0.55, size = 2, colour = GROUP_COLS[["A"]])
      }
    }

    if (use_f) p <- p + scale_colour_manual(values = GROUP_COLS)
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
