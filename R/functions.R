".cleanClass<-"<-function(x,value) {
  if (! value %in% class(x))
     class(x)<-c(value,class(x))
  x
}


.standardize_df <- function(df) {

  num <- sapply(df, is.numeric)

  df[num] <- lapply(df[num], function(x) {
    (x - mean(x, na.rm = TRUE)) / stats::sd(x, na.rm = TRUE)
  })

  df
}
