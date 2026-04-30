#' Summary for ubcBiasDetection data.frame
#'
#' This function summarize results from [ubc_lm()] of class [ubcBiasDetection]
#'
#' @param x A table of class [ubcBiasDetection]

#' @export
summary.ubcBiasDetection<-function(x,...) {

  mod<-lm(con.estimate~dis.estimate,x)
  r2<-summary(mod)$r.squared
  a<-coef(mod)[1]
  b<-coef(mod)[2]
  print(x[,c("var","adj.est","adj.se","adj.t","adj.p")])
  cat("\nTarget (winners) variables:\n\n")
  x$infl  <-  x$dis.estimate   -  x$con.estimate
  x$bias1 <-  x$dis.estimate   -  x$adj.est
  x$bias2 <-  x$con.estimate   -  x$adj.est
  df1 <- x[x$role=="W",,drop=FALSE]
  df1 <-  df1[,c("var","adj.est","infl","bias1","bias2")]
  print(df1)
  cat("\n\nControl (loosers) variables:\n\n")
  df2<-x[x$role=="L",,drop=FALSE]
  df2 <-  df2[,c("var","adj.est","infl","bias1","bias2")]
  print(df2)
  cat("\n\nOverall indices:\n\n")
  cat("Avg Inflation (S1-S2): ", mean(x$infl),"\n")
  cat("Target variables Stage 1 Inflation (S1-S2): ", mean(df1$infl),"\n")
  cat("Target variables Stage 1 bias SD: ", sd(df1$bias1),"\n")
  cat("Target variables Stage 2 bias SD: ", sd(df1$bias2),"\n")
  cat("Stages coherence: ", r2,"\n")
  cat("Stages coef (b): ", b,"\n")
  cat("Stages offest (a): ", a,"\n")
  cat("Estimation method: ", attr(x,"method"),"\n")
  res<-list(
    winners=df1,
    loosers=df2,
    a.infl=mean(x$infl),
    w.infl=mean(df1$infl),
    a.bias1.sd=sd(df1$bias1),
    a.bias2.sd=sd(df1$bias2),
    stages.r2=r2,
    stages.a=a,
    stages.b=b,
    method=attr(x,"method")
  )
  class(res)<-c("ubc.bd.summary",class(res))
  invisible(res)
}

#' Diagnostic plot for UBC adjusted estimates
#'
#' This function plots results from [ubc_lm()] of class [ubcBiasDetection]
#'
#' @param x A table of class [ubcBiasDetection]

#' @export
plot.ubcBiasDetection<-function(x,...) {

  ggplot2::ggplot(
    x,
    ggplot2::aes(
      x = estimate,
      y = c_estimate,
      color = role
    )
  ) +
    ggplot2::geom_point(size=3) +
    ggplot2::scale_color_discrete(
      name = NULL,
      labels = c(L = "Loosers", W = "Winners")
    ) +
     ggplot2::labs(
      x = "Discovery Estimates",
      y = "Stage 2 Estimates"
    ) +
    ggplot2::theme_classic(base_size = 16)
}


