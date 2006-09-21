PARmodel <- function(D, x, w = NULL, fmla = fmla, Var = c("none","delta","boot","bayes","jackknife"), CI = c("none","normal","percentile","BCa"), alpha = 0.05, B = 500)
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
    w       <- rep(1, n)
    namenp  <- colnames(x)
    if (is.null(namenp))
    stop(" colnames are missing for ", sQuote("x"))


PAR         <- rep(0,ncol(x))
VarPar      <- rep(0,ncol(x))
varsummand2 <- rep(0,ncol(x))

### Construction of C-matrix for the PAR ### 
Ci          <- bincombinations(ncol(x)-1) 
index       <- rep(0,ncol(x)-1) 
for (i in 1:length(index)) 
index[i]    <- 2^(i-1) 
ix          <- matrix(rep(sort(index,decreasing=TRUE),nrow(Ci)),ncol=length(index),byrow=TRUE) 
jplus       <- rowSums(Ci) 

    for (parti in 1:ncol(x))
    {
    x           <- cbind(x[,-1],x[,1])
    colnames(x) <- c(namesp[-1],namesp[1])
    
    ### Computation of the conditional probabilities ###
    R           <- 2^(ncol(x)-1)
    C           <- x[,-ncol(x)]
    glm.out     <- glm(formula = as.formula(fmla),"binomial",weights = w, data = data.frame(cbind(D,x)))
    preddat     <- bincombinations(ncol(x)+1)[(2^(ncol(x)+1)/2+1):(2^(ncol(x)+1)),c(1,seq((ncol(x)+1),2))]
    colnames(preddat) <- c(" ",colnames(x))
    pbed        <- matrix(predict.glm(glm.out, newdata = data.frame(preddat), type = "response"),ncol = 2)
    diffbed     <- pbed[,2]-pbed[,1]
  
    ### Construction of C-Matrix ###
    Clist       <- vector(mode = "list", length = ncol(x))
        for (m in 1:length(Clist))
        Clist[[m]]   <- x[,m]
    irg         <- interaction(Clist)
    Cs          <- model.matrix(~irg-1)
    dummy       <- bincombinations(ncol(C))
    
    
    ### Matrixmultiplikation ### 
    fak         <- gamma((ncol(x)-jplus-1)+1)*gamma(jplus+1) 
    Cit         <- t(Ci) 
    sumj        <- matrix(0,ncol=1,nrow=2^(ncol(x)-1)) 
    sumjvar     <- matrix(0,ncol=1,nrow=2^(ncol(x)-1)) 
        for (m in 1:ncol(Cit)) 
        { 
        Cij     <- rowSums((t(Cit*Ci[m,]))*ix)+1 
        sumj[m,]<- sum(fak*diffbed[Cij]) 
        sumjvar[m,] <- sum(fak*diffbed[Cij]^2) 
        } 
    sumi        <- (colSums(Cs*w)/n)[(R+1):(length(diffbed)*2)] 
    PAR[parti]  <- 1/(gamma(ncol(x)+1)*mean(D*w))*sum(sumi*sumj) 

################################################ END OF COMPUTATION OF PAR ###################################################################################
    
    if (Var=="delta" && CI=="BCa")
    stop("Computation of BCa-intervals is only possible with bootstrap or bayesian bootstrap replications!")
    if (Var=="jackknife" && CI=="BCa")
    stop("Computation of BCa-intervals is only possible with bootstrap or bayesian bootstrap replications!")

################################################ BEGIN OF COMPUTATION OF VARIANCE ESTIMATES #####################################################################

### based on delta method ###
    if (Var=="delta")
    {
    ### extraction of all main and interaction terms from glm.out$coeff ###
    betas       <- glm.out$coeff[2:(ncol(x)+1)]
    gammas      <- glm.out$coeff[(ncol(x)+2):length(glm.out$coeff)]
   
    ### Construction of design-matrix for A_mod where all coefficients from glm.out are represented in dummy-notation ###
    dummygeneral    <- cbind(rep(1,2^ncol(C)),bincombinations(ncol(x)))
    dummybeta       <- cbind(rep(1,2^ncol(C)),rep(1,2^ncol(C)),dummy)
    dummyohnebeta   <- cbind(rep(1,2^ncol(C)),rep(0,2^ncol(C)),dummy)
    colnames(dummygeneral) <- c(" ",names(betas))
    colnames(dummybeta) <- c(" ",names(betas))
    colnames(dummyohnebeta) <- c(" ",names(betas))
    
     if (length(glm.out$coeff)!=(ncol(x)+1))
     {
     helpgeneral    <- matrix(,ncol = length(gammas),nrow = nrow(dummygeneral))  
     helpbeta       <- matrix(,ncol = length(gammas),nrow = nrow(dummybeta))
     helpohnebeta   <- matrix(,ncol = length(gammas),nrow = nrow(dummyohnebeta))
     attach(data.frame(dummygeneral))
        for(gen in 1:ncol(helpgeneral))
        helpgeneral[,gen]  <- eval(parse(text = gsub(":","*",names(gammas)))[gen])
     detach(data.frame(dummygeneral))
     attach(data.frame(dummybeta))
        for(bet in 1:ncol(helpbeta))
        helpbeta[,bet]     <- eval(parse(text = gsub(":","*",names(gammas)))[bet])
     detach(data.frame(dummybeta))
     attach(data.frame(dummyohnebeta))
        for(ohn in 1:ncol(helpohnebeta)) 
        helpohnebeta[,ohn] <- eval(parse(text = gsub(":","*",names(gammas)))[ohn])
     detach(data.frame(dummyohnebeta))
     dummygeneral   <- cbind(dummygeneral,helpgeneral)
     dummybeta      <- cbind(dummybeta,helpbeta)
     dummyohnebeta  <- cbind(dummyohnebeta,helpohnebeta)
     }
    
    
    ### Calculation of A_mod of partial derivatives with respect to beta ###
    pis         <- colSums(Cs*w)/n
    Varbed      <- as.vector(pbed)*(1-as.vector(pbed))
    sumjA       <- matrix(0,ncol=length(glm.out$coeff),nrow=2^(ncol(x)-1)) 
        for (m in 1:ncol(Cit)) 
        { 
        Cij      <- rowSums((t(Cit*Ci[m,]))*ix)+1 
        sumjA[m,]<- colSums(fak*(dummybeta[Cij,]*(pbed[,2][Cij]*(1-pbed[,2][Cij]))-dummyohnebeta[Cij,]*(pbed[,1][Cij]*(1-pbed[,1][Cij]))))  
        } 
   
    Abetamod    <- -(colSums(dummygeneral*pis*Varbed))/mean(D)+colSums(sumi*sumjA)/sum(sumi*sumj)
   
    ### Computation of matrix A_mf_pi of partial derivatives w.r.t. p_i ###
    Amfpnull    <- -pbed[,1]/mean(D)
    Amfpeins    <- -pbed[,2]/mean(D)+sumj/sum(sumi*sumj)
    Amfpi       <- c(Amfpnull,Amfpeins) 
    
    ### Computation of Var(beta) ###
    Varbeta <- vcov(glm.out)*n
    
    ### Variance estimate ###
    TeilA               <- Abetamod%*%Varbeta%*%Abetamod
    TeilB               <- pis*(1-pis)*Amfpi^2
    jminus1index        <- 1-diag(length(Amfpi))
    TeilC               <- 2*pis*Amfpi*apply(pis*Amfpi*jminus1index,2,sum)
    varsummand2         <- sum(TeilB-TeilC)
    VarPar[parti]       <- PAR[parti]^2*(TeilA+sum(TeilB-TeilC))/n
    }
    } ### end of parti ###

### based on Bootstrap ###
    if (Var=="boot")
    {
    PARb            <- bootPAR(D, x, type = "boot", model = TRUE, fmla = fmla, B = B) 
    VarPar          <- apply(PARb,2,var)
    }
    
    if (Var=="bayes")
    {
    PARb            <- bootPAR(D, x, type = "bayes", model = TRUE, fmla = fmla, B = B)
    VarPar          <- apply(PARb,2,var)
    }
    
    if (Var=="jackknife")
    {
    PARb            <- bootPAR(D, x, type = "jackknife", model = TRUE, fmla = fmla)
    VarPar          <- matrix(0,ncol = ncol(x),nrow = 1) 
        for (v in 1: ncol(PARb))
        VarPar[,v]  <- (n-1)/n*sum((PARb[,v]-mean(PARb[,v]))^2)
    }
################################################ END OF COMPUTATION OF VARIANCE ESTIMATES #####################################################################
    
################################################ BEGIN OF COMPUTATION OF CONFIDENCE INTERVAL ESTIMATES ########################################################
### normal intervals ###
    if (Var!="none"  && CI=="normal")
    KI              <- matrix(c(PAR+qnorm(alpha/2)*sqrt(VarPar),PAR+qnorm(1-alpha/2)*sqrt(VarPar)),nrow = 2, byrow = TRUE)
    
    else if (Var=="boot" && CI=="normal" || Var=="jackknife" && CI=="normal")
    {
    KI              <- matrix(,ncol = ncol(x), nrow = 2)
    for (int in 1:ncol(x))
    KI[,int]        <- matrix(c(mean(PARb[,int])+qnorm(alpha/2)*sqrt(VarPar[int]),mean(PARb[,int])+qnorm(1-alpha/2)*sqrt(VarPar[int])),nrow = 2, byrow = TRUE)
    }
    
### percentile intervals
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
names(PAR)      <- namepar
names(VarPar)   <- namevar
return(list(PAR = PAR,VarPAR = VarPar))
}
else if (Var!="none" && CI!="none")
{
names(PAR)      <- namepar
names(VarPar)   <- namevar
colnames(KI)    <- nameki
rownames(KI)    <- c("CI_low","CI_up")
return(list(PAR = PAR, VarPAR = VarPar,CIPAR = KI))
}


}
