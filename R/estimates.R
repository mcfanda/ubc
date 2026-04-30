.compute_adj_se <- function(mod, newrow, dis_se, data_for_tau2,
                            vcov_type = "HC3") {

  cf <- stats::coef(mod)
  a_hat <- unname(cf[1])
  b_hat <- unname(cf[2])
  x0 <- newrow$dis.estimate[1]

  adj_est <- a_hat + b_hat * x0

  V <- tryCatch(
    sandwich::vcovHC(mod, type = vcov_type),
    error = function(e) stats::vcov(mod)
  )

  V2 <- V[1:2, 1:2, drop = FALSE]
  g <- c(1, x0)
  var_calibration <- as.numeric(t(g) %*% V2 %*% g)

  var_full <- var_calibration + (b_hat^2) * (dis_se^2)

  tau2 <- max(
    stats::var(data_for_tau2$dis.estimate, na.rm = TRUE) -
      mean(data_for_tau2$dis.SE^2, na.rm = TRUE),
    0
  )


  if (tau2 <= 0) {
    stage1_post_var <- 0
  } else {
    stage1_post_var <- (tau2 * dis_se^2) / (tau2 + dis_se^2)
  }

  var_posterior <- var_calibration + (b_hat^2) * stage1_post_var

  data.frame(
    adj.est = adj_est,
    se_calibration = sqrt(var_calibration),
    se_full = sqrt(var_full),
    se_posterior = sqrt(var_posterior),
    var_calibration = var_calibration,
    var_full = var_full,
    var_posterior = var_posterior,
    tau2 = tau2,
    stringsAsFactors = FALSE
  )
}


.loo_beta_estimates <- function(data, method, loo = TRUE, vcov_type = "HC3",
                                se_method = c("full",
                                              "posterior",
                                              "calibration",
                                              "none")) {

  se_method <- match.arg(se_method)

  if (!inherits(data, "ubcBiasDetection")) {
    stop("UBC method requires data of class `ubcBiasDetection`")
  }

  if (!requireNamespace("sandwich", quietly = TRUE)) {
    stop("Package 'sandwich' is required")
  }

  if (nrow(data) < 3) {
    warning("Too few coefficients to estimate unbiased results")
    return(data.frame(
      adj.est = rep(NA_real_, nrow(data)),
      adj.se = rep(NA_real_, nrow(data)),
      adj.var = rep(NA_real_, nrow(data)),
      var.calibration = rep(NA_real_, nrow(data)),
      var.stage1 = rep(NA_real_, nrow(data)),
      tau2 = rep(NA_real_, nrow(data)),
      adj.t = rep(NA_real_, nrow(data)),
      adj.p = rep(NA_real_, nrow(data)),
      stringsAsFactors = FALSE
    ))
  }

  if (loo && nrow(data) < 10) {
    loo <- FALSE
  }

  fit_model <- function(d) {
    if (toupper(method) == "WLS") {
      stats::lm(
        con.estimate ~ dis.estimate,
        data = d,
        weights = 1 / (d$con.SE^2)
      )
    } else {
      stats::lm(
        con.estimate ~ dis.estimate,
        data = d
      )
    }
  }

  pick_se <- function(se_row) {

    if (se_method == "none") {
      adj_se <- NA_real_
      adj_var <- NA_real_
      var_stage1 <- NA_real_
    } else if (se_method == "calibration") {
      adj_se <- se_row$se_calibration
      adj_var <- se_row$var_calibration
      var_stage1 <- 0
    } else if (se_method == "full") {
      adj_se <- se_row$se_full
      adj_var <- se_row$var_full
      var_stage1 <- se_row$var_full - se_row$var_calibration
    } else if (se_method == "posterior") {
      adj_se <- se_row$se_posterior
      adj_var <- se_row$var_posterior
      var_stage1 <- se_row$var_posterior - se_row$var_calibration
    } else {
      stop("Unknown se_method")
    }

    data.frame(
      adj.est = se_row$adj.est,
      adj.se = adj_se,
      adj.var = adj_var,
      var.calibration = se_row$var_calibration,
      var.stage1 = var_stage1,
      tau2 = se_row$tau2,
      stringsAsFactors = FALSE
    )
  }

  out_list <- vector("list", nrow(data))

  if (loo) {

    for (i in seq_len(nrow(data))) {

      d_fit <- data[-i, , drop = FALSE]
      mod <- fit_model(d_fit)

      se_row <- .compute_adj_se(
        mod = mod,
        newrow = data[i, , drop = FALSE],
        dis_se = data$dis.SE[i],
        data_for_tau2 = d_fit,
        vcov_type = vcov_type
      )

      out_list[[i]] <- pick_se(se_row)
    }

  } else {

    mod <- fit_model(data)

    for (i in seq_len(nrow(data))) {

      se_row <- .compute_adj_se(
        mod = mod,
        newrow = data[i, , drop = FALSE],
        dis_se = data$dis.SE[i],
        data_for_tau2 = data,
        vcov_type = vcov_type
      )

      out_list[[i]] <- pick_se(se_row)
    }
  }

  out <- do.call(rbind, out_list)
  rownames(out) <- NULL

  out$adj.t <- out$adj.est / out$adj.se
  out$adj.p <- 2 * stats::pnorm(-abs(out$adj.t))

  out
}

### helper


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


