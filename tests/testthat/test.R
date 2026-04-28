simulate_correlated_y <- function(n, r, seed = NULL, offset=0) {

  if (!is.null(seed)) {
    set.seed(seed)
  }

  if (any(abs(r) > 1)) {
    stop("All target correlations must be between -1 and 1.")
  }

  k <- length(r)

  y <- stats::rnorm(n)
  X <- matrix(NA_real_, nrow = n, ncol = k)

  for (j in seq_len(k)) {
    e <- stats::rnorm(n)
    X[, j] <- r[j] * y + sqrt(1 - r[j]^2) * e
  }
  varnames<-c("y",names(r))
  if (is.null(names(r))) varnames<- c("y", paste0("x", seq_len(k)))
  out <- data.frame(y = y, X, check.names = FALSE)
  names(out) <- varnames

  out[] <- lapply(out, function(z) as.numeric(scale(z, center = TRUE, scale = TRUE))+offset)

  return(out)
}

r<-c(rep(.4,20),rep(0,25))
names(r)<-paste0("x",1:45)
data<-simulate_correlated_y(100,r)
disc<-ubc::lm_candidates(y~.,data)
r2<-r[disc$var]
data2<-simulate_correlated_y(100,r2)
res<-ubc::ubc_lm(y~.,disc,data2)
res
summary(res)
