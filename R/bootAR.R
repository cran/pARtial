bootAR <- function(D, x, C = NULL,  model = NULL, fmla = fmla,  type = c("boot","bayes","jackknife"), B = 500) 
{
type    <- match.arg(type)
n       <- length(D)

if (type=="boot")
{
    bs      <- rmultinom(B, n, rep(1, n)/n)
    ARb     <- rep(0, B)
        for (b in 1:B)
        {
            if (is.null(model))
            ARb[b]  <- AR(D = D, x = x, C = C, w = bs[,b])$ARisk
            else if (!is.null(model))
            ARb[b]  <- AR(D, x, C, model = TRUE, fmla = fmla, w = bs[,b])$ARisk
        }
}

else if (type=="bayes")
{
    ARb  <- rep(0, B)
    
    agam    <- rep(1,n)
    l       <- length(agam)
    gam     <- matrix(rgamma(l * B, agam), ncol = l, byrow = TRUE)
    sm      <- gam %*% rep(1, l)
    bs      <- gam/as.vector(sm)
    
        for ( b in 1: B)
        {
            if (is.null(model))
            ARb[b] <- AR(D, x, C, w = round(bs[b,]*n))$ARisk
            else if (!is.null(model))
            ARb[b] <- AR(D, x, C, model = TRUE, fmla = fmla, w = as.vector(round(bs[b,]*n)))$ARisk
        }
}

else if (type=="jackknife")
{   
    index   <- rep(1,n)
    ARb     <- rep(0,n)
        for(j in 1:length(index))
        {
        bs      <- index
        bs[j]   <- 0
            if (is.null(model))
            ARb[j]  <- AR(D, x, C, w = bs)$ARisk
            else if (!is.null(model))
            ARb[j]  <- AR(D, x, C, model = TRUE, fmla = fmla, w = bs)$ARisk
        }
}   

return(ARb)

}
