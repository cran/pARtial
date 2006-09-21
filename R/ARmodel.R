ARmodel <- function(D , x , C = NULL, w = NULL, fmla, Var = c("none","delta","boot","bayes","jackknife"), CI = c("none","normal","logit","percentile","BCa"), B = 500, alpha = 0.05) 
{

if (is.vector(x))
x       <- matrix(x,ncol=1)
n       <- length(D)
namenx  <- colnames(x)
namenC  <- colnames(C)

if (!is.null(C))
    {
    if (is.vector(C))
    C   <- matrix(C,ncol = 1)
    if (!is.matrix(C) || nrow(C) != n) 
        stop(paste(sQuote("C"), " is not a matrix with ", n, " rows"))
    }
    if (is.matrix(D))
        stop(paste(sQuote("D"), " is not a vector"))
    if (nrow(x) != n)
        stop(paste(sQuote("x"), " is not of length ", n))
    if (is.null(w)) 
    w   <- rep(1, n)
    else if (!is.null(w))
    w   <- w  
    
### unadjusted attributable risk ###

if (is.null(C))
{   glm.out     <- glm(formula = as.formula(fmla), "binomial", data = data.frame(x), weights = w) 
    preddat     <- bincombinations(ncol(x)+1)[(2^(ncol(x)+1)/2+1):2^(ncol(x)+1),]
    colnames(preddat) <- c(" ",namenx)
    pbed        <- predict.glm(glm.out, newdata = data.frame(preddat), type = "response")
    Attrib      = 1-(pbed[1]/mean(D*w))
}

### adjusted attributable risk ###
else if (!is.null(C))
{
    glm.out     <- glm(formula = as.formula(fmla), "binomial", data = data.frame(cbind(x,C)), weights = w)
    preddat     <- bincombinations(ncol(x)+ncol(C)+1)[(2^(ncol(x)+ncol(C)+1)/2+1):2^(ncol(x)+ncol(C)+1),]
    colnames(preddat) <- c(" ",namenx,namenC)
    pbed        <- predict.glm(glm.out, newdata = data.frame(preddat), type = "response")
    pbed2       <- matrix(pbed*(1-pbed),ncol=1)
    pbedohneE   <- pbed[1:(2^ncol(C))]
    pbedohneE2  <- pbed2[1:(2^ncol(C))]
    pbedmitE    <- pbed[(2^ncol(C)+1):nrow(preddat)]
    pbedmitE2   <- pbed2[(2^ncol(C)+1):nrow(preddat)]
        
### estimation of P(C) and P(D) ###    
   
    Clist       <- vector(mode = "list", length = ncol(C))
    for (m in 1:length(Clist))
    Clist[[m]]  <- C[,m]
    irg         <- interaction(Clist)
    Cs          <- model.matrix(~irg-1)    
    pC          <- (colSums(Cs*w)/n)
    sortindex   <- sort(names(pC),index.return=TRUE)
    pC          <- pC[sortindex$ix]
    pC          <- ifelse(pC==0,pC+1e-30,pC)
    pD          <- mean(D*w)
    pD          <- ifelse(pD==0,pD+1e-30,pD)
    
### Attributable Risk ###
    Z           <- sum(pbedohneE*pC)
    Attrib      <- 1-(Z/pD)

}

############################################### End of point estimation ################################################
if (Var=="delta")
{
    comb        <- bincombinations(ncol(x))
    ### Covariance matrix of coefficients from glm ###
    Vb          <- vcov(glm.out)

############################################## Variance of unadjusted AR ###############################################
    if (is.null(C)) 
    {    
    ### Calculation of variance of multinomial random vector sigma ###
    Cplist      <- vector(mode = "list", length = ncol(x)+1)
    Cplist[[1]] <- D
    for(li in 2:length(Cplist))
    Cplist[[li]]<- x[,li-1]
    irgp        <- interaction(Cplist)
    Cp          <- model.matrix(~irgp-1)
    pij         <- colSums(Cp*w)/n
    sortind     <- sort(names(pij),index.return=TRUE)
    pij         <- pij[sortind$ix]
    px          <- pij[1:(length(pij)/2)]+pij[(length(pij)/2)+1:length(pij)][1:(length(pij)/2)]
    sigma       <- (diag(pij)-pij%*%t(pij))/n
    sigma1      <- (diag(px)-px%*%t(px))/n
    ### Computation of Hessian matrix ###
     betas      <- glm.out$coeff[(ncol(x)+2):length(glm.out$coeff)] ### extraction of all interaction terms from glm.out$coeff if ncol(x)!=1 ###
     helpmatrix <- cbind(rep(1,nrow(comb)),comb)
     colnames(helpmatrix) <- c(" ",namenx)
     attach(data.frame(helpmatrix))
     if(length(glm.out$coeff)==(ncol(x)+1))
        H       <- t(rbind(helpmatrix*(-n*pbed),helpmatrix*(n-n*pbed)))
     else if (length(glm.out$coeff)!=(ncol(x)+1))
        {
        helpmatrix2      <- matrix(,ncol = length(betas),nrow = nrow(comb))
        for (pm in 1:ncol(helpmatrix2))
        helpmatrix2[,pm] <- eval(parse(text = gsub(":","*",names(betas)))[pm])
        helpmatrix2      <- cbind(helpmatrix,helpmatrix2)
        H                <- t(rbind(helpmatrix2*(-n*pbed),helpmatrix2*(n-n*pbed)))
        }
     detach(data.frame(helpmatrix))
    ### Computation of variance covariance matrix of p and b ###
    Chelp       <- Vb%*%H%*%sigma
    Cvar        <- Chelp[,1:(ncol(Chelp)/2)]+Chelp[,((ncol(Chelp)/2)+1):ncol(Chelp)]
   
    ### Calculation of partial derivatives with respect to b (JA) and p(BA) ###
    pbed2       <- pbed*(1-pbed)
    JAalpha     <- pbed2[1]/pbed[1]-sum(pbed2*px)/sum(pbed*px)
        if(ncol(x)==1)
        JAbeta           <- -(pbed2[2]*px[2])/sum(pbed*px)
        else if (ncol(x)!=1 && length(glm.out$coeff)==(ncol(x)+1))
        JAbeta           <- -colSums(comb*pbed2*px)/sum(pbed*px)
        else if (ncol(x)!=1 && length(glm.out$coeff)!=(ncol(x)+1))
        JAbeta           <- -colSums(helpmatrix2[,-1]*pbed2*px)/sum(pbed*px) 
    JA          <- matrix(c(JAalpha,JAbeta))
    BA          <- -pbed/sum(pbed*px)
    
    ### Computation of variance estimate ###
    Varar       <- (1-Attrib)^2*((t(JA)%*%Vb%*%JA)+(t(BA)%*%sigma1%*%BA)+(2*t(BA)%*%t(Cvar)%*%JA))
    }

############################################################# Variance estimate for adjusted AR #######################################################

    else if (!is.null(C)) 
    {
    comb2       <- bincombinations(ncol(x)+ncol(C))
    comb3       <- bincombinations(ncol(C))
   
    ### Calculation of variance of multinomial random vector sigma ###
    pijk        <- as.vector(t(matrix(prop.table(ftable(data.frame(D,x,C))),ncol=2)))
    pCx         <- pijk[1:(length(pijk)/2)]+pijk[(length(pijk)/2)+1:length(pijk)][1:(length(pijk)/2)]
    sigma       <- (diag(pijk)-pijk%*%t(pijk))/n
    sigma1      <- (diag(pCx)-pCx%*%t(pCx))/n
   
    ### Computation of Hessian matrix ###
    gammas      <- glm.out$coeff[(ncol(x)+ncol(C)+2):length(glm.out$coeff)] ### extraction of all interaction terms from glm.out$coeff ###
    helpmatrix  <- cbind(rep(1,nrow(comb2)),comb2)  
    colnames(helpmatrix) <- c(" ",namenx,namenC)
    attach(data.frame(helpmatrix))
        if(length(glm.out$coeff)==(ncol(x)+ncol(C)+1))
        H           <- t(rbind(helpmatrix*(-n*pbed),helpmatrix*(n-n*pbed)))
        else if (length(glm.out$coeff)!=(ncol(x)+ncol(C)+1))
        {
        helpmatrix2 <- matrix(,ncol = length(gammas),nrow = nrow(comb2))
        for (pm in 1:ncol(helpmatrix2))
        helpmatrix2[,pm] <- eval(parse(text=gsub(":","*",names(gammas)))[pm])
        helpmatrix2 <- cbind(helpmatrix,helpmatrix2)
        H           <- t(rbind(helpmatrix2*(-n*pbed),helpmatrix2*(n-n*pbed)))
        }
    
    ### Computation of variance covariance matrix of p and b ###
    Cov         <- Vb%*%H%*%sigma
    Cvar        <- Cov[,1:(ncol(Cov)/2)]+Cov[,((ncol(Cov)/2)+1):ncol(Cov)]

    ### Calculation of JA and BA ###
    if(length(glm.out$coeff)==(ncol(x)+ncol(C)+1))
    {
    JA1         <- colSums(cbind(rep(1,nrow(comb3)),matrix(0,ncol=ncol(x),nrow=nrow(comb3)),comb3)*pbedohneE2*as.vector(pC))/sum(pbedohneE*as.vector(pC))
    JA2         <- colSums(cbind(rep(1,nrow(comb2)),comb2)*as.vector(pbed2)*as.vector(pCx))/sum(pbed*as.vector(pCx))
    JA          <- matrix(JA1-JA2,nrow=1)
    }
    
    else if(length(glm.out$coeff)!=(ncol(x)+ncol(C)+1))
    {
    helpmatrix3 <- matrix(,ncol = length(gammas), nrow = nrow(helpmatrix))
    for(xi in 1:ncol(helpmatrix3))
    helpmatrix3[,xi] <- eval(parse(text=gsub(":","*",names(gammas)))[xi])
    helpmatrix3 <- cbind(helpmatrix,helpmatrix3)[1:(ncol(C)*2),]
    JA1         <- colSums(helpmatrix3*pbedohneE2*as.vector(pC))/sum(pbedohneE*as.vector(pC))
    JA2         <- colSums(helpmatrix2*as.vector(pbed2)*as.vector(pCx))/sum(pbed*as.vector(pCx))
    JA          <- matrix(JA1-JA2,nrow=1)
    }
    BA          <- matrix(pbedohneE/sum(pbedohneE*as.vector(pC))-pbed/sum(pbed*as.vector(pCx)),nrow=1)

    ### Calculation of variance estimates ###
    Sep         <- sqrt((1-Attrib)^2*(JA%*%Vb%*%t(JA)))
    Varar       <- (1-Attrib)^2*((JA%*%Vb%*%t(JA))+(BA%*%sigma1%*%t(BA))+(2*BA%*%t(Cvar)%*%t(JA)))
    }

}
    else if(Var=="boot")
    {
    ARb     <- bootAR(D, x, C , type = "boot", model = TRUE, fmla = fmla)
    Varar   <- var(ARb)
    }
    else if(Var=="bayes")
    {
    ARb     <- bootAR(D, x, C , type="bayes", model = TRUE, fmla = fmla)
    Varar   <- var(ARb)
    }
    else if(Var=="jackknife")
    {
    ARb     <- bootAR(D, x, C , type="jackknife", model = TRUE, fmla = fmla)
    Varar   <- (n-1)/n*sum((ARb-mean(ARb))^2)
    } 


################################################ BEGIN OF COMPUTATION OF CONFIDENCE INTERVAL ESTIMATES ########################################################

### normal intervals ###
if (Var=="delta"  && CI=="normal")
KI          <- c(Attrib+qnorm(alpha/2)*sqrt(Varar),Attrib+qnorm(1-alpha/2)*sqrt(Varar))
else if (Var=="boot" && CI=="normal" ||Var=="bayes" && CI=="normal"||Var=="jackknife" && CI=="normal")
KI          <- c(mean(ARb)+qnorm(alpha/2)*sqrt(Varar),mean(ARb)+qnorm(1-alpha/2)*sqrt(Varar))

### logit-transformed normal intervals ###
else if (Attrib==0 && !is.null(CI)&& CI=="logit")
stop(paste("AR=",Attrib,":","The logit-transformation is not feasible if AR=",Attrib,"!"))
else if (Var=="delta" && CI=="logit") 
KI          <- c(1/(1+((1-Attrib)/Attrib)*exp((-qnorm(alpha/2)*sqrt(Varar))/(Attrib*(1-Attrib)))),1/(1+((1-Attrib)/Attrib)*exp((qnorm(alpha/2)*sqrt(Varar))/(Attrib*(1-Attrib)))))
else if (Var=="boot" && CI=="logit" ||Var=="bayes" && CI=="logit"||Var=="jackknife" && CI=="logit")
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
else if (Var=="boot" && CI=="percentile" || Var=="bayes" && CI=="percentile")
{
    ARsort  <- sort(ARb)
    KI      <- c(ARsort[B*(alpha/2)],ARsort[B*(1-alpha/2)])
}

### BCa intervals ###
else if (Var=="boot" && CI=="BCa" ||Var=="bayes" && CI=="BCa")
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

if(Var=="none")
return(list(ARisk = Attrib))

else if (Var!="none" && CI=="none")
return(list(ARisk = Attrib,VarARisk = Varar, SEARisk = sqrt(Varar)))

else if (Var!="none" && CI!="none")
return(list(ARisk = Attrib, VarARisk = Varar,CIARisk = KI))

}
