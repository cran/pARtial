bootPAR <- function(D, x, model = NULL, fmla = fmla, type = c("boot","bayes","jackknife"), B = 500) 
{
type    <- match.arg(type)
namesp  <- names(data.frame(x))
n       <- length(D)
if (type=="boot")
{
    bs          <- rmultinom(B, n, rep(1, n)/n)
    PARb        <- matrix(0,ncol = ncol(x),nrow = B)
    for (b in 1:B)
    {
        if (is.null(model))
        PARb[b,]     <- PartialAR(D, x, w = as.vector(bs[,b]))
        else if (!is.null(model))
        PARb[b,]     <- PartialAR(D, x, model = TRUE, fmla = fmla, w = as.vector(bs[,b]))
    }        
    
}

else if (type=="bayes")
{
    agam        <- rep(1,n)
    l           <- length(agam)
    gam         <- matrix(rgamma(l * B, agam), ncol = l, byrow = TRUE)
    sm          <- gam %*% rep(1, l)
    bs          <- gam/as.vector(sm)#

    PARb        <- matrix(0,ncol = ncol(x), nrow = B)
    for (b in 1:B)
    {
        if (is.null(model))
        PARb[b,]    <- PartialAR(D,x,w = as.vector(round(bs[b,]*n)))
        else if (!is.null(model))
        PARb[b,]    <- Parmod(D, x, w = as.vector(round(bs[b,]*n)),fmla = fmla)
    }    
}

else if (type=="jackknife")
{   
    PARb        <- matrix(0,ncol = ncol(x),nrow = n)
    for(j in 1:n)
    {
        if (is.null(model))
        PARb[j,]    <- PartialAR(D[-j], x[-j,])
        else if (!is.null(model))
        PARb[j,]    <- PartialAR(D[-j], x[-j,], model = TRUE, fmla = fmla)
    }    
}    

namepar <- rep(0,ncol(x))
for(nam in 1:ncol(x))
    {
        if(is.null(namesp))
        namepar[nam]    <- paste("Rep_PAR(E",c(1:ncol(x))[nam],")",sep="")    
        else if (!is.null(namesp))
        namepar[nam]    <- paste("Rep_PAR(",namesp[nam],")",sep="")
    }
colnames(PARb) <- namepar
return(PARb)

}
    
   
     
