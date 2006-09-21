PARmodelfree = function(D, x, w = NULL, model = NULL, fmla, Var = c("none","delta","boot","bayes","jackknife"), CI = c("none","normal","percentile","BCa"), alpha = 0.05, B = 500)
{
Var <- match.arg(Var)
CI  <- match.arg(CI)

n   <- length(D)
namesp <- names(data.frame(x))
if (!is.matrix(x))
stop(sQuote("x"), " must be a matrix ")
if (!is.vector(D))
stop(sQuote("D"), " must be a vector")
if (is.null(w)) 
w  <- rep(1, n)

PAR         <- rep(0,ncol(x))
VarPar      <- rep(0,ncol(x))
varsummand2 <- rep(0,ncol(x))

### Construction of C-matrix for the PAR ### 
CiLminus1   <- bincombinations(ncol(x)-1) 
index       <- rep(0,ncol(x)-1) 
for (i in 1:length(index)) 
index[i]    <- 2^(i-1) 
ix          <- matrix(rep(sort(index,decreasing=TRUE),nrow(CiLminus1)),ncol=length(index),byrow=TRUE) 
jplus       <- rowSums(CiLminus1) 

    for (parti in 1:ncol(x)) 
    { 
    x           <- cbind(x[,-1],x[,1]) 
    ### computation of conditional probabilities ### 
    R           <- 2^(ncol(x)-1) 
    Clist       <- vector(mode = "list", length = ncol(x)) 
        for (m in 1:length(Clist)) 
        Clist[[m]]   <- x[,m] 
    irg         <- interaction(Clist) 
    Cs          <- model.matrix(~irg-1) 
    pbed        <- colSums(D*Cs*w)/colSums(Cs*w) 
    pbed        <- ifelse(pbed=="NaN",0,pbed) 
    pbed        <- ifelse(pbed=="Inf",0,pbed) 
    pbedA       <- pbed[1:R] 
    pbedB       <- pbed[(R+1):length(pbed)] 
    diffbed     <- pbedB-pbedA 

    ### Multiplication of matrices ### 
    fak         <- gamma((ncol(x)-jplus-1)+1)*gamma(jplus+1) 
    Cit         <- t(CiLminus1) 
    sumj        <- matrix(0,ncol=1,nrow=2^(ncol(x)-1)) 
    sumjvar     <- matrix(0,ncol=1,nrow=2^(ncol(x)-1)) 
        for (m in 1:ncol(Cit)) 
        { 
        Cij         <- rowSums((t(Cit*CiLminus1[m,]))*ix)+1    ### Cij ist Index, welche bedingten Ws. in jeder j-ten Schicht genommen werden müssen!###
        sumj[m,]    <- sum(fak*diffbed[Cij]) 
        sumjvar[m,] <- sum(fak*diffbed[Cij]^2) 
        } 

    sumi            <- (colSums(Cs*w)/n)[(R+1):length(pbed)] 
    
    PAR[parti]      <- 1/(gamma(ncol(x)+1)*mean(D*w))*sum(sumi*sumj) 
    PAR[parti]      <- ifelse(PAR[parti]=="NaN",0,PAR[parti])
    PAR[parti]      <- ifelse(PAR[parti]=="Inf",0,PAR[parti])

################################################ END OF COMPUTATION OF PAR ###################################################################################
    
    if (Var=="delta" && CI=="BCa")
    stop("Computation of BCa-intervals is only possible with bootstrap or bayesian bootstrap replications!")
    if (Var=="jackknife" && CI=="BCa")
    stop("Computation of BCa-intervals is only possible with bootstrap or bayesian bootstrap replications!")

################################################ BEGIN OF COMPUTATION OF VARIANCE ESTIMATES #####################################################################

### based on delta method ###
    if (Var=="delta")
    {
    Pi          <- colSums(Cs*w)/n
    ### Computation of matrix A_mf of partial derivatives w.r.t. m_i ###
    Amfauf      <- matrix(0,ncol=1,nrow=2^(ncol(x)-1))
    sumjv       <- matrix(0,ncol=1,nrow=2^(ncol(x)-1))
        for (Ci in 1:nrow(CiLminus1))
        {
            for (mv in 1:ncol(Cit))
            {
            Cijv        <- rowSums((t(Cit*CiLminus1[mv,]))*ix)+1
            Cijv        <- ifelse(Cijv==CiLminus1,1,0)
            sumjv[mv,]  <- sum(sumi[mv]*sum(fak*Cijv))
            }
        Amfauf[Ci,]  <- sum(sumjv)
        }
    Amf         <- (-Pi/mean(D))+(c(-Amfauf,Amfauf)/sum(sumi*sumj))
 
    ### Variance of  mu_i and p_i ###
    Varmui      <- (pbed*(1-pbed))/Pi
    Varmui      <- ifelse(Varmui=="NaN",0,Varmui)
    Varmui      <- ifelse(Varmui=="Inf",0,Varmui)
    Varpi       <- diag(Pi)-Pi%*%t(Pi)
    
   
    ### Computation of matrix A_mf_pi of partial derivatives w.r.t. p_i ###
    Amfpnull    <- -pbedA/mean(D)
    Amfpeins    <- -pbedB/mean(D)+sumj/sum(sumi*sumj)
    Amfpi       <- c(Amfpnull,Amfpeins) 
    
    ### Variance estimate ###
    TeilA               <- sum(Varmui*Amf^2)
    TeilB               <- Pi*(1-Pi)*Amfpi^2
    jminus1index        <- 1-diag(length(Amfpi))
    TeilC               <- 2*Pi*Amfpi*apply(Pi*Amfpi*jminus1index,2,sum)
    VarPar[parti]       <- PAR[parti]^2*(TeilA+sum(TeilB-TeilC))/n
    
    }
    } ### end of parti-loop ###


### based on Bootstrap ###
    if (Var=="boot")
    {
    PARb            <- bootPAR(D, x, type = "boot", B = B) 
    VarPar          <- apply(PARb,2,var)
    }
    
    if (Var=="bayes")
    {
    PARb            <- bootPAR(D, x, type = "bayes", B = B)
    VarPar          <- apply(PARb,2,var)
    }
    
    if (Var=="jackknife")
    {
    PARb            <- bootPAR(D, x, type = "jackknife")
    VarPar          <- matrix(0,ncol = ncol(x),nrow = 1) 
        for (v in 1: ncol(PARb))
        VarPar[,v]  <- (n-1)/n*sum((PARb[,v]-mean(PARb[,v]))^2)
    }
################################################ END OF COMPUTATION OF VARIANCE ESTIMATES #####################################################################
    
################################################ BEGIN OF COMPUTATION OF CONFIDENCE INTERVAL ESTIMATES ########################################################

### normal intervals ###
    if (Var=="delta"  && CI=="normal")
    KI              <- matrix(c(PAR+qnorm(alpha/2)*sqrt(VarPar),PAR+qnorm(1-alpha/2)*sqrt(VarPar)),nrow = 2, byrow = TRUE)
    
    else if (Var=="boot" && CI=="normal" || Var=="jackknife" && CI=="normal")
    {
    KI              <- matrix(,ncol = ncol(x), nrow = 2)
    for (int in 1:ncol(x))
    KI[,int]        <- matrix(c(mean(PARb[,int])+qnorm(alpha/2)*sqrt(VarPar[int]),mean(PARb[,int])+qnorm(1-alpha/2)*sqrt(VarPar[int])),nrow = 2, byrow = TRUE)
    }

### percentile intervals ###
    else if (Var=="boot" && CI=="percentile" ||Var=="bayes" && CI=="percentile")
    {
    KI              <- matrix(,ncol = ncol(x), nrow = 2)
    PARsort         <- apply(PARb,2,sort)
        for(int in 1:ncol(x))
        KI[,int]    <- c(PARsort[,int][B*(alpha/2)],PARsort[,int][B*(1-alpha/2)])
    }

### BCa-intervals ###
    else if (Var=="boot" && CI=="BCa" ||Var=="bayes" && CI=="BCa")
    {
    KI              <- matrix(,ncol = ncol(x), nrow = 2)     
    resultjack      <- bootPAR(D , x, type = "jackknife")
    PARsort         <- apply(PARb,2,sort)
    
        for (int in 1:ncol(x))
        {
        counts      <- length(which(PARsort[,int] <= PAR[int]))
        counts      <- ifelse(counts==0,counts+1e-7,counts)
        counts      <- ifelse(counts==B,counts-1e-4,counts)
        z0          <- qnorm(counts/B)    
        a           <- sum((mean(resultjack[,int])-resultjack[,int])^3)/(6*(sum((mean(resultjack[,int])-resultjack[,int])^2))^1.5)
        a           <- ifelse(a=="NaN",0,a)
        alpha1      <- pnorm(z0+((z0+qnorm(alpha/2))/(1-a*(z0+qnorm(alpha/2)))))
        alpha2      <- pnorm(z0+((z0+qnorm(1-alpha/2))/(1-a*(z0+qnorm(1-alpha/2)))))   
         if (B*alpha1<1 && B*alpha2>1)
         KI[,int]   <- c(PARsort[,int][1],PARsort[,int][B*alpha2])
         else if (B*alpha2<1 && B*alpha1>1)
         KI[,int]   <- c(PARsort[,int][B*alpha1],PARsort[,int][1])
         else if (B*alpha1<1 && B*alpha2<1)
         KI[,int]   <- rep(PARsort[,int][1],2)
         else
         KI[,int]   <- c(PARsort[,int][B*alpha1],PARsort[,int][B*alpha2])
         }    
    }
################################################ END OF COMPUTATION OF CONFIDENCE INTERVAL ESTIMATES ###############################################################

####################################### Preparation of print out of results ################################################################
namenhelp       <- 1:ncol(x)
namepar         <- vector(length=ncol(x))
namevar         <- vector(length=ncol(x))
nameki          <- vector(length=ncol(x))
    for(nam in 1:ncol(x))
    {
        if(is.null(namesp))
        {
        namepar[nam]    <- paste("PAR(E",namenhelp[nam],")",sep="")    
        namevar[nam]    <- paste("Var(PAR(E",namenhelp[nam],"))_",Var,sep="")  
        nameki[nam]     <- paste("CI(PAR(E",namenhelp[nam],"))_",CI,sep="")  
        }
        else if (!is.null(namesp))
        {
        namepar[nam]    <- paste("PAR(",namesp[nam],")",sep="")
        namevar[nam]    <- paste("Var(PAR(",namesp[nam],"))",sep="")  
        nameki[nam]     <- paste("CI(PAR(",namesp[nam],"))",sep="")  
        }
    }
  
if(Var=="none")
{
names(PAR) <- namepar
return(PAR = PAR)
}
else if (Var!="none" && CI=="none")
{
names(PAR) <- namepar
names(VarPar)   <- namevar
return(list(PAR = PAR,VarPAR = VarPar))
}
else if (Var!="none" && CI!="none")
{
names(PAR) <- namepar
names(VarPar)   <- namevar
colnames(KI)    <- nameki
rownames(KI)    <- c("CI_low","CI_up")
return(list(PAR = PAR, VarPAR = VarPar,CIPAR = KI))
}

}
