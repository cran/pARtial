ARmodelfree <- function(D, x, C = NULL, w = NULL, Var = c("none","delta","boot","bayes","jackknife"), CI = c("none","normal","logit","percentile","BCa"), alpha = 0.05, B = 500) 
{

if (!is.vector(x) && ncol(x)!=1)
x       <- as.vector(Strata(x))
else if (!is.vector(x) && ncol(x)==1)
x       <- as.vector(x)
n       <- length(D)


if (!is.null(C))
    {
    if (is.vector(C))
    C   <- matrix(C,ncol = 1)
    if (!is.matrix(C) || nrow(C) != n) 
        stop(paste(sQuote("C"), " is not a matrix with ", n, " rows"))
    }
    if (is.matrix(D))
        stop(paste(sQuote("D"), " is not a vector"))
    if (length(x) != n)
        stop(paste(sQuote("x"), " is not of length ", n))
    if (is.null(w)) 
    w   <- rep(1, n)


if (is.null(C))
Attrib  <- 1-(((sum(w*D)-sum(x*D*w))/(n-sum(w*x)))/(sum(w*D)/n))

   
else if (!is.null(C))
{
    Clist   <- vector(mode = "list", length = ncol(C))
        for (m in 1:length(Clist))
        Clist[[m]]   <- C[,m]
        
    irg     <- interaction(Clist)
    Cs      <- model.matrix(~irg-1)    
    
    wneu    <- w / sum(w)
    frac1   <- (colSums(D*Cs*wneu) - colSums(D*Cs*x*wneu))
    frac2   <- (colSums(Cs*wneu) - colSums(x * Cs * wneu))
    frac2   <- ifelse(frac2==0,frac2+1e-30,frac2)
    frac    <- frac1/frac2
    pD      <- sum(D*wneu)
    pD      <- ifelse(pD==0,pD+1e-30,pD)
    if(pD==0)
    Attrib  <- 0
    else if (pD!=0)
    Attrib  <- 1 - sum(frac * colSums(Cs * wneu)) / pD
}
################################################ END OF COMPUTATION OF AR ###################################################################################


################################################ BEGIN OF COMPUTATION OF VARIANCE ESTIMATES #####################################################################
if (Var =="delta")
{   if(is.null(C))
    {
    w       <- w/sum(w)
    x11     <- sum(D*x*w)
    x12     <- sum(x*w)-x11
    x21     <- sum(D*w)-sum(D*x*w)
    x22     <- 1-(x11+x12+x21)
    Varar   <- x21*(x11*x22*(1-x21)+(x21^2)*x12)/(n*((x11+x21)^3)*((x21+x22)^3))
    }
    
    if(!is.null(C))
    {
    phi     <- -Attrib+1
    p11     <- colSums(x*D*Cs*w)
    p12     <- colSums(x*Cs*w)-p11
    p21     <- colSums(D*Cs*w)-p11
    p22     <- colSums(Cs*w)-(p11+p12+p21)
    p11k    <- colSums(x*D*Cs)/colSums(Cs)
    p12k    <- (colSums(x*Cs)-colSums(x*D*Cs))/colSums(Cs)
    p21k    <- (colSums(D*Cs)-colSums(x*D*Cs))/colSums(Cs)
    p22k    <- (colSums(Cs)-(colSums(x*D*Cs)+(colSums(x*Cs)-colSums(x*D*Cs))+(colSums(D*Cs)-colSums(x*D*Cs))))/colSums(Cs)
    p.1.    <- sum(p11+p21)
    p2.k    <- p21k+p22k
    p.1k    <- p11k+p21k
    p..k    <- colSums(Cs*w)
    n..k    <- colSums(Cs)
            
    A1      <- (phi/p.1.)^2
    A1      <- ifelse(A1=="Inf",0,A1)
    A1      <- ifelse(A1=="NaN",0,A1)
    A2      <- p..k^2*p11k/n..k
    A2      <- ifelse(A2=="Inf",0,A2)
    A2      <- ifelse(A2=="NaN",0,A2)
    AA      <- A1*sum(A2)
            
    B1      <- (p..k/p.1.)^2
    B1      <- ifelse(B1=="NaN",0,B1)
    B1      <- ifelse(B1=="Inf",0,B1)
    B2      <- (p22k/p2.k^2-phi)^2*p21k/n..k
    B2      <- ifelse(B2=="NaN",0,B2)
    B2      <- ifelse(B2=="Inf",0,B2)
    BB      <- sum(B1*B2)
            
    C1      <- (1/p.1.)^2
    C1      <- ifelse(C1=="NaN",0,C1)
    C1      <- ifelse(C1=="Inf",0,C1)
    C2      <- p..k^2*p21k^2*p22k/(n..k*p2.k^4)
    C2      <- ifelse(C2=="NaN",0,C2)
    C2      <- ifelse(C2=="Inf",0,C2)
    CC      <- C1*sum(C2)
            
    D1      <- p..k^2*p.1k^2/n..k
    D1      <- ifelse(D1=="NaN",0,D1)
    D1      <- ifelse(D1=="Inf",0,D1)
    D2      <- sum(D1)/p.1.^2
    D2      <- ifelse(D2=="NaN",0,D2)
    D2      <- ifelse(D2=="Inf",0,D2)
    DD      <- phi^2*D2
            
    Varar   <- abs(AA+BB+CC-DD)
    }
}
if(Var=="boot")
{
    ARb     <- bootAR(D, x, C ,  type="boot")
    Varar   <- var(ARb)
}
if(Var=="bayes")
{
    ARb     <- bootAR(D, x, C , type="bayes")
    Varar   <- var(ARb)
}
if(Var=="jackknife")
{
    ARb     <- bootAR(D, x, C , type="jackknife")
    Varar   <- (n-1)/n*sum((ARb-mean(ARb))^2)
}

################################################ BEGIN OF COMPUTATION OF CONFIDENCE INTERVAL ESTIMATES ########################################################
    if (Var=="delta" && CI=="BCa")
    stop("Computation of BCa-intervals is only possible with bootstrap or bayesian bootstrap replications!")
    if (Var=="jackknife" && CI=="BCa")
    stop("Computation of BCa-intervals is only possible with bootstrap or bayesian bootstrap replications!")
### normal intervals ###
if (Var=="delta"  && CI=="normal")
KI          <- c(Attrib+qnorm(alpha/2)*sqrt(Varar),Attrib+qnorm(1-alpha/2)*sqrt(Varar))
if (Var=="boot" && CI=="normal" ||Var=="bayes" && CI=="normal"||Var=="jackknife" && CI=="normal")
KI          <- c(mean(ARb)+qnorm(alpha/2)*sqrt(Varar),mean(ARb)+qnorm(1-alpha/2)*sqrt(Varar))

### logit-transformed normal intervals ###
if (Attrib==0 && CI=="logit")
stop(paste("AR=",Attrib,":","The logit-transformation is not feasible if AR=",Attrib,"!"))
if (Var=="delta" &&  CI=="logit") 
KI          <- c(1/(1+((1-Attrib)/Attrib)*exp((-qnorm(alpha/2)*sqrt(Varar))/(Attrib*(1-Attrib)))),1/(1+((1-Attrib)/Attrib)*exp((qnorm(alpha/2)*sqrt(Varar))/(Attrib*(1-Attrib)))))
if (Var=="boot" && CI=="logit" ||Var=="bayes" && CI=="logit"||Var=="jackknife" && CI=="logit")
    {
    ptrunc      <- length(which(ARb <= 0))/length(ARb)
            if(ptrunc > 0.95)
            stop("proportion of negative replications is over 95%")
    ARbtrunc    <- ARb[ARb > 0]
    ARtlogit    <- log(ARbtrunc/(1-ARbtrunc))
    Etlogit     <- mean(ARtlogit) 
    if (Var!="jackknife")
    Vartlogit   <- var(ARtlogit)
    else if (Var=="jackknife")
    Vartlogit   <- (n-1)/n*sum((ARtlogit-mean(ARtlogit))^2)         
                if(ptrunc == 0)
                {
                ARb     <- ifelse(ARb==1,1-1e-5,ARb)
                ARlogit <- log(ARb/(1-ARb))
                if (Var!="jackknife")
                Varlogit<- var(ARlogit)
                if (Var=="jackknife")
                Varlogit<- (n-1)/n*sum((ARlogit-mean(ARlogit))^2)
                Cil     <- c(mean(ARlogit)+qnorm(alpha/2)*sqrt(Varlogit),mean(ARlogit)+qnorm(1-alpha/2)*sqrt(Varlogit))
                KI      <- exp(Cil)/(1+exp(Cil))
                }
                
                else if (ptrunc != 0 && ptrunc <= 0.95)
                {
                Alpha       <- qnorm(ptrunc)
                lambda      <- dnorm(Alpha)/(1-pnorm(Alpha))
                delta       <- lambda*(lambda-Alpha)
                sigmalogit  <- Vartlogit/(1-delta)
                mulogit     <- Etlogit-sqrt(sigmalogit)*lambda
                Cil         <- c(mulogit+qnorm(alpha/2)*sqrt(sigmalogit),mulogit-qnorm(alpha/2)*sqrt(sigmalogit))            
                KI          <- exp(Cil)/(1+exp(Cil))
                }

    }
    
### percentile intervals ###
if (Var=="boot" && CI=="percentile" || Var=="bayes" && CI=="percentile")
{
    ARsort  <- sort(ARb)
    KI      <- c(ARsort[B*(alpha/2)],ARsort[B*(1-alpha/2)])
}

### BCa intervals ###
if (Var=="boot" && CI=="BCa" ||Var=="bayes" && CI=="BCa")
{
    ARsort  <- sort(ARb)
    counts  <- 0  
            for (ar in 1:length(ARb))
            {
            if (ARb[ar]<Attrib)
            counts  <- counts+1
            }
    counts  <- ifelse(counts==0,counts+1e-7,counts)
    counts  <- ifelse(counts==B,counts-1e-4,counts)
    z0      <- qnorm(counts/B)    
    resultjack  <- bootAR(D , x, C, type = "jackknife") 
    a       <- sum((mean(resultjack)-resultjack)^3)/(6*(sum((mean(resultjack)-resultjack)^2))^1.5)
    a       <- ifelse(a=="NaN",0,a)
    alpha1  <- pnorm(z0+((z0+qnorm(alpha/2))/(1-a*(z0+qnorm(alpha/2)))))
    alpha2  <- pnorm(z0+((z0+qnorm(1-alpha/2))/(1-a*(z0+qnorm(1-alpha/2)))))
         if (B*alpha1<1 && B*alpha2>1)
         KI      <- c(ARsort[1],ARsort[B*alpha2])
         else if (B*alpha2<1 && B*alpha1>1)
         KI      <- c(ARsort[B*alpha1],ARsort[1])
         else if (B*alpha1<1 && B*alpha2<1)
         KI      <- rep(ARsort[1],2)
         else
         KI      <- c(ARsort[B*alpha1],ARsort[B*alpha2])
}


################################################ END OF COMPUTATION OF CONFIDENCE INTERVAL ESTIMATES ###############################################################

if(Var=="none" && CI=="none")
return(list(ARisk = Attrib))

if (Var!="none"&& CI=="none")
return(list(ARisk = Attrib,VarARisk = Varar))

if (Var!="none" && CI!="none")
return(list(ARisk = Attrib, VarARisk = Varar,CIARisk = KI))
 
}
