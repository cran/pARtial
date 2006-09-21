boot <- function(D, x, C = NULL, param = c("AR","PAR"), model = NULL, fmla = fmla, type = c("boot","bayes","jackknife"), B = 500) 
{
type <- match.arg(type)
param <- match.arg(param)

if (param=="PAR" && is.null(C)==FALSE)
stop("Matrix with variables for adjustment are only valid for computation of adjusted AR!")


if (param=="AR")
Rep <- bootAR(D = D, x = x, C = C, model = model, fmla = fmla, type = type, B = B)

else if (param=="PAR")
Rep <- bootPAR(D = D, x = x, model = model, fmla = fmla, type = type, B = B)

return(Rep)
    
}
