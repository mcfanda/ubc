#' ubcBiasDetection objects
#'
#' Objects of class `ubcBiasDetection` are returned by [ubc_lm()].
#'
#' They contain the estimated parameters and related statistics.
#'
#' @section Structure:
#' Objects of class `ubcBiasDetection` are lists with the following elements:
#'
#'  A `ubcBiasDetection` object must contain at least the following columns:
#'
#' \describe{
#'   \item{var}{the name of the candidate variable}
#'   \item{c_estimate}{Discovery (Stage 1) estimates}
#'   \item{estimate}{Stage 2 estimates}
#'   \item{adj.est}{UBC adjusted estimates}
#'   \item{role}{`W` for winner variables and `L` for control loosers variables}
#' }
#' @name ubcBiasDetection
NULL
