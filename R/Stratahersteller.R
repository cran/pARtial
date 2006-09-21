Strata <- function(x) 
{

    Neu     <- matrix(,ncol=1,nrow=dim(x)[1])
    Neu     <- cbind(Neu,x)
    
    for (k in 2:ncol(Neu))
    {
    Neu[,1][Neu[,k]==1]<-1
    }
    
    xneu    <- matrix(Neu[,1],ncol=1)
    xneu[is.na(xneu)]<-0

return(xneu)



}
