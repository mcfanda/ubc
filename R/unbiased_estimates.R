#' Find the variables candidated for selections (winners)
#'
#' This function estimates a general model parameters for a set
#' of variables and identifies the plausible candidates for selection (winners) and
#' a set of control variables (loosers) based on significance level
#'
#' @param formula A formula of the form y~. or y~x+. in which the dot `.` is replaced
#'  by each input variable. A GLM is estimated for each `varnames` variable by replace it in the formula
#' @param data A data.frame containing the variables to analyze
#' @param varnames a vector of variables names. If NULL, all variables in the dataframe are tested, excluding the dependent variable
#' @param index the index of the coefficient in the model results to analyze. Default is 2, indicating the first effect in a GLM (1 would be the intercept)
#' @param sig.level The critical alpha level to be used in the selection of the variable, default .10
#' @param prop.loser Proportion of control variables (loosers) to select with respect of the selected ones (winners). Default .33 indicates selecting 1 looser for 3 winners
#' @return A table with the estimates, SE, t-test and p-values of the selected variables
#' @export
lm_candidates <- function(formula, data, varnames = NULL, index = 2,
                           sig.level = .10, prop.loser = .33) {

  ests <- ubc_lm_estimates(formula, data, varnames = varnames, index = index)

  if (is.null(ests) || nrow(ests) == 0) {
    warning("No sensible variables found for candidate selection")
    out <- .empty_lm_candidates()
    .cleanClass(out) <- "ubcResults"
    return(out)
  }

  winners <- ests[ests$p < sig.level, , drop = FALSE]
  .cleanClass(winners) <- "ubcResults"
  W <- nrow(winners)
  E <- nrow(ests)

  if (W == 0) {
    warning("No candidate found as winners")
    out <- ests[0, , drop = FALSE]
    out$role <- character(0)
    .cleanClass(out) <- "ubcResults"
    return(out)
  }

  winners$role <- "W"

  if (W == E) {
    warning("All candidates are selected as winners")
    return(winners)
  }

  L <- round(prop.loser * W)
  loosers <- ests[ests$p >= sig.level, , drop = FALSE]
  L1 <- nrow(loosers)

  if (L1 == 0)
    return(winners)

  if (L1 < L)
    L <- L1

  candidates <- winners

  if (L > 0) {
    loosers <- loosers[sample.int(nrow(loosers), L), , drop = FALSE]
    loosers$role <- "L"
    candidates <- rbind(winners, loosers)
  }

  .cleanClass(candidates) <- "ubcResults"
  candidates
}


#' Estimate coefficients of a GLM
#'
#' This function estimates standard coefficients of a general model based useful for the UCB method.
#' Coefficients are computed considering applying the model template to all input variables.
#'
#' @param formula A formula of the form y~. or y~x+. in which the dot `.` is replaced
#'  by each input variable. A GLM is estimated for each `varnames` variable by replace it in the formula
#' @param data A data.frame
#' @param index the index of the coefficient in the model results to analyze. Default is 2, indicating the first effect in a GLM (1 would be the intercept)
#' @return A data.frame with the unbiased estimates, SE, t-test and p-values of the selected variables
#' @export
ubc_lm_estimates <- function(formula, data, varnames = NULL, index = 2) {

  if (!inherits(formula, "formula"))
    stop("Please provide a model formula of class `formula`")

  if (!is.data.frame(data))
    stop("`data` must be a data.frame")

  lhs <- all.vars(formula)[1]

  if (!lhs %in% names(data))
    stop("Dependent variable not found in `data`")

  if (is.null(varnames))
    varnames <- setdiff(names(data), lhs)

  varnames <- unique(varnames)
  varnames <- varnames[varnames %in% names(data)]

  if (length(varnames) == 0) {
    warning("No candidate variables available in `data`")
    return(.empty_lm_estimates())
  }

  forms <- .make_formulas(formula, varnames)

  mods <- lapply(seq_along(forms), function(i) {
    .safe_lm_estimate(forms[[i]], data = data, varname = varnames[i], index = index)
  })

  ok <- !vapply(mods, is.null, logical(1))

  if (!any(ok)) {
    warning("No sensible variables could be estimated")
    return(.empty_lm_estimates())
  }

  if (any(!ok)) {
    warning(sprintf(
      "%s variable(s) could not be estimated and were dropped: %s",
      sum(!ok),
      paste(varnames[!ok], collapse = ", ")
    ))
  }

  results <- as.data.frame(do.call(rbind, mods[ok]), stringsAsFactors = FALSE)
  row.names(results) <- NULL
  results
}


#' Estimate unbiased coefficients of a GLM
#'
#' This function estimates unbiased coefficients of a general model based on the UCB method.
#' Coefficients are computed considering both the discovery data and the stage 2 data.
#'
#' @param formula A formula of the form y~. or y~x+. in which the dot `.` is replaced
#'  by each input variable. A GLM is estimated for each `varnames` variable by replace it in the formula
#' @param candidates Discovery (stage 1) data with candidates variables (winners and controlling loosers).
#' @param stage2 A data.frame with stage 2 (confirmatory) data
#' @param index the index of the coefficient in the model results to analyze. Default is 2, indicating the first effect in a GLM (1 would be the intercept)
#' @param method `"WLS"` (default) weights the estimation with 1/SE^2, where SE is the Stage 2 estimates standard errors.
#'                 `"OLS"` does not apply any weights.
#' @param loo apply the leave-one-out algorithm to estimate the unbiased coefficients
#' @return A data.frame with the unbiased estimates, SE, t-test and p-values of the selected variables
#' @export
ubc_lm <- function(formula, candidates, stage2, varnames = NULL, index = 2, method = "WLS", loo=TRUE) {

  if (!"ubcResults" %in% class(candidates))
    stop("Please provide discovery phase results of class ubcResults")

  if (is.null(candidates) || nrow(candidates) == 0) {
    warning("There are no candidate variables in input")
    return(NULL)
  }

  if (is.null(varnames))
    varnames <- candidates$var

  stage2 <- .standardize_df(stage2)
  confirmation <- ubc_lm_estimates(formula, stage2, varnames = varnames, index = index)

  if (is.null(confirmation) || nrow(confirmation) == 0) {
    warning("No sensible stage-2 estimates could be computed")
    return(NULL)
  }

  names(confirmation) <- paste0("c_", names(confirmation))
  names(confirmation)[1] <- "var"

  cand <- candidates
  cand_role <- cand$role
  cand$role <- NULL

  data <- merge(confirmation, cand, by = "var", all.x = TRUE, sort = FALSE)

  if (nrow(data) < 3) {
    warning("Too few candidates to estimate unbiased estimates")
    return(NULL)
  }

  idx <- match(data$var, candidates$var)
  data$role <- candidates$role[idx]

  .cleanClass(data) <- "ubcBiasDetection"
  data$adj.est <- .loo_beta_estimates(data, method = method,loo=loo)
  attr(data, "method") <- method
  attr(data, "loo") <- loo

  data
}


### unexported functions

.loo_beta_estimates <- function(data, method,loo=TRUE) {

  if (!"ubcBiasDetection" %in% class(data))
    stop("UBC method requires data of class `ubcBiasDetection`")

  if (nrow(data) < 3) {
    warning("Too few rows for unbiased estimation")
    return(rep(NA_real_, nrow(data)))
  }

  if (loo && nrow(data) < 5) {
    warning("Too few rows for leave-one-out estimation")
    return(rep(NA_real_, nrow(data)))
  }

  myweights <- 1 / (data$c_SE^2)

  if (loo) {
  coefs <- lapply(seq_len(nrow(data)), function(i) {

    .data <- data[-i, , drop = FALSE]

    .weights <- NULL
    if (toupper(method) == "WLS")
      .weights <- myweights[-i]

    out <- tryCatch({
      mod <- stats::lm(c_estimate ~ estimate, data = .data, weights = .weights)
      stats::predict(mod, newdata = data[i, , drop = FALSE])[1]
    }, error = function(e) {
      NA_real_
    }, warning = function(w) {
      invokeRestart("muffleWarning")
    })

    as.numeric(out)
  })

  results<-unlist(coefs)
  } else {
    .weights <- NULL
    if (toupper(method) == "WLS")
      .weights <- myweights

    mod <- stats::lm(c_estimate ~ estimate, data = data, weights = .weights)
    coefs<-stats::predict(mod)
    results<-unlist(coefs)

}

  return(results)

}

.safe_lm_estimate <- function(formula, data, varname, index) {

  out <- tryCatch({

    mod <- stats::lm(formula, data = data)
    cf  <- summary(mod)$coefficients

    if (is.null(dim(cf)))
      return(NULL)

    if (nrow(cf) < index)
      return(NULL)

    vals <- cf[index, ]

    data.frame(
      var = varname,
      estimate = unname(vals[1]),
      SE = unname(vals[2]),
      t = unname(vals[3]),
      p = unname(vals[4]),
      stringsAsFactors = FALSE
    )

  }, error = function(e) {
    NULL
  }, warning = function(w) {
    invokeRestart("muffleWarning")
  })

  out
}


.empty_lm_estimates <- function() {
  data.frame(
    var = character(0),
    estimate = numeric(0),
    SE = numeric(0),
    t = numeric(0),
    p = numeric(0),
    stringsAsFactors = FALSE
  )
}


.empty_lm_candidates <- function() {
  data.frame(
    var = character(0),
    estimate = numeric(0),
    SE = numeric(0),
    t = numeric(0),
    p = numeric(0),
    role = character(0),
    stringsAsFactors = FALSE
  )
}


.make_formulas <- function(formula, vars) {
  lapply(vars, function(v) {
    f_new <- do.call(
      substitute,
      list(formula, list(`.` = as.name(v)))
    )
    as.formula(f_new, env = environment(formula))
  })
}
