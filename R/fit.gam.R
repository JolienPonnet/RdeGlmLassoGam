##' Function used to determine initial GAM estimate
##' @param Ry the response vector.
##' @param RB the model matrix.
##' @param RrS penalty matrix.
##' @param Rfamily  a character string indicating the family. This can be "poisson" or "binomial".
fit.gam.sp1 <- function (Ry, RB, RrS, Rfamily)
.Call("fit_gam_sp_cpp", Ry, RB, RrS, Rfamily, PACKAGE = "RdeGlmLassoGam")

