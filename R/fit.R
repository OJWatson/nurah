#' Fit DAG-based Bayesian Hierarchical Model
#'
#' This function constructs and fits a Bayesian hierarchical linear model (via **brms** and Stan)
#' to estimate "indirect mortality" given a user-specified causal DAG and observed data. It uses
#' the DAG to determine which predictors (direct causes of the outcome) to include, and adds
#' specified spatial hierarchy levels as random intercept effects. The modeling approach
#' follows principles from McElreath's *Statistical Rethinking* book, using DAGs to inform model
#' structure and employing weakly-informative priors and partial pooling for robust inference.
#'
#' @param data A data frame of observed data. Must contain the `outcome` variable and all predictor and grouping variables.
#' @param dag A `nurah_dag` object defining the causal structure. This DAG is used to extract the direct causes of `outcome` as fixed-effect predictors.
#' @param outcome Name (string) of the dependent variable in `data` to model.
#' @param spatial_levels Character vector of spatial grouping variables for hierarchical random effects. For example, `c("region", "district")` indicates a multilevel structure with districts nested within regions (first element is the highest level).
#' @param priors Optional brms prior specification. If `NULL`, defaults to weakly informative priors: `Normal(0, 0.5)` for fixed-effect coefficients and `HalfNormal(0, 1)` for standard deviation parameters. (The half-normal is implemented by using a `Normal(0,1)` prior on `sd` class parameters, which are inherently non-negative).
#' @param chains Number of MCMC chains to run (defaults to 4).
#' @param cores Number of CPU cores to use for parallel sampling (defaults to `getOption("mc.cores", 1)`; set higher for parallel chains).
#' @param iter Number of iterations per chain (including warmup; default 2000).
#' @param warmup Number of warmup (burn-in) iterations per chain (default is half of `iter` if not specified).
#' @param seed Random seed for reproducibility (optional).
#' @param control A list of control parameters passed to Stan (e.g., `list(adapt_delta = 0.95)` to help with convergence).
#' @param formula_override Optional formula object to use as the model formula, overriding the default formula construction from the DAG. If provided, this formula is used as-is (the user should include any desired random effects terms in it).
#' @param family The likelihood family for the outcome, specified as a brms family object. Defaults to `gaussian()` (normal) for continuous outcomes; can be changed to, e.g., `binomial()` or `poisson()` depending on the `outcome`.
#' @param ... Additional arguments passed to `brms::brm()`, such as `sample_prior`, `control` settings, etc.
#'
#' @details
#' **Model Construction:** By default, the fixed-effect formula is derived from the provided DAG.
#' The DAG's structure is used to include only the direct causes (parents) of the outcome node as
#' predictors in a linear model for the outcome. This ensures a principled, causally-informed model
#' that avoids adding extraneous variables or adjusting for mediators/colliders without causal justification. All such predictor variables must exist as columns in `data`.
#'
#' **Hierarchical Structure:** The `spatial_levels` argument allows adding random intercepts for
#' specified grouping factors to model multi-level variation (e.g., variation between regions and
#' districts). If multiple levels are provided, they are assumed to be nested in the order given
#' (e.g., `region/district`). Including hierarchical random effects induces partial pooling across
#' these groups, which guards against overfitting by shrinking extreme group-level estimates toward the overall mean. This reflects the idea that groups with less data get stronger shrinkage, a key benefit of multilevel models described in *Statistical Rethinking*.
#'
#' **Priors:** In line with *Statistical Rethinking* recommendations, the function uses conservative
#' default priors to regularize estimates. Population-level coefficients get a Normal(0, 0.5) prior,
#' a mildly informative prior assuming predictors are on roughly standard scales. Standard deviations
#' of group-level effects (and the residual error) get a half-normal prior (implemented as Normal(0,1)
#' on the standard deviation). These priors provide gentle skepticism toward very large effects or variances, while still being weak enough not to overwhelm the data.
#'
#' **Missing Data:** Rows with missing values in the outcome or any of the predictor or grouping variables used are dropped before fitting, since `brms` (and the underlying Stan sampler) cannot handle `NA` values.
#'
#' **Model Fitting:** The model is fit using Hamiltonian Monte Carlo via `brms::brm()`. By default, a
#' Gaussian (normal) likelihood is used for a continuous outcome; this can be changed by specifying the
#' `family` argument (for example, use `binomial()` for binary mortality outcomes). You can pass further
#' arguments to `brm()` via `...` for additional control (e.g., setting `iter`, `warmup`, or `adapt_delta`
#' for more complex models). The result is a `brmsfit` object containing the posterior samples and model details,
#' which can be examined with standard methods (`summary()`, `plot()`, posterior predictive checks, etc.).
#'
#' @return A `brmsfit` object containing the fitted Bayesian model (result of the `brm()` call).
#'
#' @examples
#' \dontrun{
#' # Define a causal DAG: suppose 'mortality' is directly caused by 'nutrition' and 'access_to_care'
#' dag <- define_dag(
#'   mortality ~ nutrition + access_to_care,
#'   nutrition ~ socioeconomic_status,
#'   access_to_care ~ socioeconomic_status + region
#' )
#'
#' # Simulate some data consistent with the DAG structure
#' set.seed(42)
#' N <- 100
#' region <- rep(1:5, each = 20)
#' district <- rep(1:20, each = 5)
#' socioeconomic_status <- rnorm(N)
#' nutrition <- 0.5 * scale(socioeconomic_status)[,1] + rnorm(N, 0, 0.1)
#' region_eff_access <- rnorm(5, 0, 0.5)       # region effect on access_to_care
#' access_to_care <- 0.5 * scale(socioeconomic_status)[,1] + region_eff_access[region] + rnorm(N, 0, 0.5)
#' region_eff_mort <- rnorm(5, 0, 1.0)       # region effect on mortality
#' district_eff_mort <- rnorm(20, 0, 0.2)      # district effect on mortality
#' mortality <- -0.3 * nutrition - 0.5 * access_to_care + region_eff_mort[region] +
#'              district_eff_mort[district] + rnorm(N, 0, 0.5)
#' data <- data.frame(region, district, socioeconomic_status, nutrition, access_to_care, mortality)
#'
#' # Fit the Bayesian model using the DAG for predictors and region/district as random effects
#' fit <- fit_dag(data = data, dag = dag, outcome = "mortality",
#'                                    spatial_levels = c("region", "district"),
#'                                    iter = 1000, chains = 2)
#' summary(fit)
#' }
#'
#' @export
fit_dag <- function(data, dag, outcome,
                                        spatial_levels = NULL,
                                        priors = NULL,
                                        chains = 4,
                                        cores = getOption("mc.cores", 1),
                                        iter = 2000,
                                        warmup = floor(iter / 2),
                                        seed = NULL,
                                        control = list(),
                                        formula_override = NULL,
                                        family = gaussian(),
                                        ...) {
  # Basic input validations
  if (!is.character(outcome) || length(outcome) != 1) {
    stop("`outcome` must be a single string (name of the outcome variable).")
  }
  if (!(outcome %in% names(data))) {
    stop("Outcome variable '", outcome, "' not found in `data`.")
  }

  # Construct fixed-effects formula from DAG (unless overridden)
  if (is.null(formula_override)) {
    if (!inherits(dag, "nurah_dag")) {
      stop("`dag` must be a 'nurah_dag' object created with define_dag().")
    }
    # Extract direct causes (parents) of the outcome from the DAG
    direct_causes <- NULL
    if ("dagitty" %in% class(dag) || !is.null(attr(dag, "dag"))) {
      # If dag is or contains a dagitty object
      if (!requireNamespace("dagitty", quietly = TRUE)) {
        stop("Package 'dagitty' is required to extract DAG information.")
      }
      # Underlying dagitty object may be the object itself or stored in an attribute
      dag_obj <- if ("dagitty" %in% class(dag)) dag else attr(dag, "dag")
      direct_causes <- tryCatch(
        dagitty::parents(dag_obj, outcome),
        error = function(e) NULL
      )
    } else if (!is.null(dag$parents)) {
      # If DAG stores parent sets in a list or environment
      direct_causes <- dag$parents[[outcome]]
    } else if (!is.null(dag$relationships)) {
      # If DAG stores formula relationships in a list
      rel <- dag$relationships[[outcome]]
      if (!is.null(rel)) {
        if (inherits(rel, "formula")) {
          # Parse formula to get right-hand side variable names
          direct_causes <- all.vars(rhs <- rel[[3]])
        } else if (is.character(rel)) {
          direct_causes <- rel
        }
      }
    }
    if (is.null(direct_causes)) {
      # If outcome not found in DAG or has no parents, use intercept-only
      warning("No parents of '", outcome, "' defined in DAG; using an intercept-only model.")
      direct_causes <- character(0)
    }
    # Build the fixed effects formula string
    if (length(direct_causes) == 0) {
      fixed_formula_str <- paste(outcome, "~ 1")
    } else {
      fixed_formula_str <- paste(outcome, "~", paste(direct_causes, collapse = " + "))
    }
    formula_fixed <- stats::as.formula(fixed_formula_str)
  } else {
    if (!inherits(formula_override, "formula")) {
      stop("`formula_override` must be a formula object (e.g., outcome ~ predictors) if provided.")
    }
    formula_fixed <- formula_override
  }

  # Construct random-effects formula for spatial hierarchy if specified
  random_formula_str <- NULL
  if (!is.null(spatial_levels) && length(spatial_levels) > 0) {
    # Ensure grouping variables exist in data
    missing_groups <- setdiff(spatial_levels, names(data))
    if (length(missing_groups) > 0) {
      stop("Grouping variable(s) not found in data: ", paste(missing_groups, collapse = ", "))
    }
    if (length(spatial_levels) > 1) {
      # Assume nested hierarchy in the given order
      random_formula_str <- paste0("(1|", paste(spatial_levels, collapse = "/"), ")")
    } else {
      random_formula_str <- paste0("(1|", spatial_levels, ")")
    }
  }

  # Combine fixed and random parts into the final formula
  if (!is.null(random_formula_str)) {
    final_formula <- stats::as.formula(paste(deparse(formula_fixed), "+", random_formula_str))
  } else {
    final_formula <- formula_fixed
  }

  # Remove rows with missing data in outcome or any required predictors/grouping variables
  vars_needed <- c(outcome)
  if (is.null(formula_override)) {
    vars_needed <- c(vars_needed, direct_causes)
  } else {
    # If formula was overridden, include all vars on its right-hand side
    vars_needed <- c(vars_needed, all.vars(formula_fixed[[3]]))
  }
  if (!is.null(spatial_levels)) {
    vars_needed <- c(vars_needed, spatial_levels)
  }
  vars_needed <- unique(vars_needed)
  data_complete <- data[stats::complete.cases(data[, vars_needed, drop = FALSE]), ]
  if (nrow(data_complete) < nrow(data)) {
    message("Dropped ", nrow(data) - nrow(data_complete), " rows due to missing values in outcome or predictors.")
  }

  # Set default priors if none provided
  if (is.null(priors)) {
    priors <- c(
      brms::prior("normal(0, 0.5)", class = "b"),        # Coefficients ~ Normal(0, 0.5)
      brms::prior("normal(0, 1)", class = "sd"),         # Random-effect SDs ~ Half-Normal(0, 1)
      brms::prior("normal(0, 1)", class = "sigma")       # Residual SD ~ Half-Normal(0, 1) for Gaussian
    )
  }

  # Fit the Bayesian model using brms (HMC via Stan)
  fit <- brms::brm(
    formula = final_formula,
    data    = data_complete,
    family  = family,
    prior   = priors,
    chains  = chains,
    cores   = cores,
    iter    = iter,
    warmup  = warmup,
    seed    = seed,
    control = control,
    ...
  )

  return(fit)
}
