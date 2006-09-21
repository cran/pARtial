AR <- function(D, x, C = NULL, model = NULL, fmla, w = NULL, Var = c("none","delta","boot","bayes","jackknife"), CI = c("none","normal","logit","percentile","BCa"), alpha = 0.05, B = 500) 
{
Var <- match.arg(Var)
CI  <- match.arg(CI)
if(any(is.na(c(D,x,C)))==TRUE)
stop("Data set contains missing values! Remove missing values with na.omit() first!")

    if (Var=="delta" && CI=="BCa" || Var=="delta" && CI=="percentile")
    stop("Computation of percentile or BCa-intervals is only possible with bootstrap or bayesian bootstrap replications!")
    if (Var=="jackknife" && CI=="BCa" || Var=="jackknife" && CI=="percentile")
    stop("Computation of percentile or BCa-intervals is only possible with bootstrap or bayesian bootstrap replications!")
if(is.null(model))
Result  <- ARmodelfree(D = D, x = x, C = C, w = w, Var = Var, CI = CI, alpha = alpha, B = B)
else if (!is.null(model))
Result  <- ARmodel(D = D, x = x, C = C, w = w, fmla = fmla, Var = Var, CI = CI, alpha = alpha, B = B)
return(Result)

} 
