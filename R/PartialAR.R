PartialAR <- function(D, x, w = NULL, model = NULL, fmla, Var = c("none","delta","boot","bayes","jackknife"), CI = c("none","normal","percentile","BCa"), alpha = 0.05, B = 500)
{
Var <- match.arg(Var)
CI  <- match.arg(CI)
if(any(is.na(c(D,x)))==TRUE)
stop("Data set contains missing values! Remove missing values with na.omit() before proceeding!")

    if (Var=="delta" && CI=="BCa")
    stop("Computation of BCa-intervals is only possible with bootstrap or bayesian bootstrap replications!")
    if (Var=="jackknife" && CI=="BCa")
    stop("Computation of BCa-intervals is only possible with bootstrap or bayesian bootstrap replications!")
if(is.null(model))
Result  <- PARmodelfree(D = D, x = x, w = w, Var = Var, CI = CI, alpha = alpha, B = B)
else if (!is.null(model))
Result  <- PARmodel(D = D, x = x, w = w, fmla = fmla, Var = Var, CI = CI, alpha = alpha, B = B)
return(Result)



}
