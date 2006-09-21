\name{AR}
\alias{AR}
\title{(Adjusted) Attributable Risk Estimates, Variances and Confidence Intervals}
\description{
AR derives crude, adjusted, crude joint or adjusted joint attributable risk estimates for one or multiple
exposure factors of primary interest adjusted for secondary confounding variables together with different
variance estimates and confidence intervals.
}
\usage{
AR( D, x, C = NULL, model = NULL, fmla, w = NULL, 
    Var = c("none","delta","boot","bayes","jackknife"), 
    CI = c("none","normal","logit","percentile","BCa"), alpha = 0.05, B = 500) 
}

\arguments{
  \item{D}{a vector containing a dichotomous indicator variable for the disease status.}
  \item{x}{a matrix containing a dichotomous indicator variable for the exposure status. If \code{ncol(x)=1}, the crude or adjusted
            attributable risk for the risk factor in \code{x} is computed. If \code{ncol(x)>1}, the joint attributable risk of the multiple
            risk factors in \code{x} is returned.}
  \item{C}{a matrix containing one or multiple confounding variables for adjusting the attributable risk. Every column of \code{C} must be dichotomous. 
            If \code{C=NULL}, the crude attributable risk for the exposure in \code{x} is computed.}
  \item{w}{a weight vector which is used to define the resampling technique. If a nonparametric or bayesian bootstrap
            or the jackknife is used, \code{w} can be ignored by the user as it is regulated by
            the input parameter \code{Var}. If else the user wants to use a different
            resampling method, \code{w} can individually be changed.    }
  \item{model}{if \code{model=TRUE}, the attributable risk is computed by use of coefficients from a logistic regression model. If
                \code{model=NULL}, the attributable risk is computed with probabilities directly estimated from the contingency tables of the data set.}
  \item{fmla}{if \code{model=TRUE}, \code{fmla} defines the desired form of the logistic regression model and is an obligatory parameter.}
  \item{Var}{a character string indicating the method of variance estimation: \code{Var="delta"} indicates
            a variance estimate derived over the delta method, \code{Var="boot"} means application of a nonparametric bootstrap, \code{Var="bayes"} indicates the Bayesian
            Bootstrap and \code{Var="jackknife"} the Jackknife. If the default \code{Var="none"} is selected, only the point estimate of the attributable risk is returned.
            \code{Var="none"} is the default!}
  \item{CI}{a character string indicating the method of confidence interval estimation: if \code{CI="normal"} a confidence
            interval constructed by using percentiles from a standard normal distribution is computed. If \code{CI="logit"} a logit-transformation of
            the attributable risk is used. If the logit-transformation is used together with variance estimation based on resampling methods (bootstrap or jackknife)
            moments of a truncated normal distribution are used for the construction of the confidence interval if the empirical distribution of the attributable risk
            contains negative values. \code{"percentile"} and \code{CI="BCa"} yields confidence intervals based on the simple percentile method and the BCa method,
            respectively (only possible when \code{Var="boot"} or \code{Var="bayes"}). \code{CI="none"} is the default!}
  \item{alpha}{the probability of error for the estimation of confidence intervals, yielding a $1-\alpha$ confidence level.}
  \item{B}{number of replications for resampling methods.}
}

\details{Point estimates of the crude attributable risk are the same wether the model--free or the model--based approach is applied.
For adjusted attributable risk estimates, the model--free approach yields estimates based on the case--load weighting approach, where the attributable risk
is written as a weighted-sum over all adjustment levels of the confounding variables in \code{C}. If the model--based approach is applied based on
a main effects model (for instance 
\code{fmla = "D~x1+c1+c2"}) the point estimate is based on the Mantel-Haenszel approach. This means that an adjusted Odds Ratio
is plugged into the formula of the attributable risk. By using a fully saturated model for estimation 
(thus \code{fmla="D~(x1+c1+c2)^3"}), the point estimate
again equals the model--free estimate, as the interaction structure within the data is totally considered. If adjusted attributable risks are computed for
large data sets containing no sparse strata, the model--free approach is recommended. But if many of the strata defined by the adjustment levels of the confounding
variables contain only few observations, the model--based approach should be used. The main benefit of the model--based approach is its flexibility, as a logistic regression 
model perfectly mapping the interaction structure of the dataset can be determined over the parameter \code{fmla}.

All supplied variants of variance estimation yield asymptotic estimates. The variance estimate based on the delta method is computationally least expensive. It is
based on a expansion of the attributable risk about its mean by taking a one step Taylor
approximation. 
The nonparametric bootstrap (\code{Var="boot"}) is based on Efron's bootstrap. Samples are taken out of the underlying data set
and the attributable risk is repeatedly computed for each of those samples. Rubin's bayesian bootstrap (\code{Var="bayes"}) is a weighted method of the bootstrap.
A weight vector is generated from a Dirichlet distribution. For each sampled weight vector, the attributable risk is computed with weighted observations.
The jackknife can be seen as an approximation to Efron's nonparametric bootstrap. Replications of the attributable risk are computed by leaving out every
single observation once at a time. If the data set contains far more observations as the number of replications normally computed for a bayesian or nonparametric
bootstrap, the computation of a bootstrap is recommended. For the number of replications, \code{B=500} is usually assumed to suffice. The number of replications should
be increased if BCa-intervals are computed. BCa-intervals thus are computationally expensive, but yield very stable results in many situations. Attention: Results of 
simulation studies showed that the application of the simple percentile method can yield unsatisfactory results when applied to the adjusted AR. 
}
  
\value{
  If attributable risk estimates are computed a single value for the point estimate is returned.
  If \code{Var!="none"} a list containing the point estimate together with its corresponding variance estimate is returned.
  If \code{CI!="none"} a list containing the point estimate, the variance estimate and the corresponding confidence interval is returned. 
}

\note{The variables \code{D}, \code{x} and \code{C} have to be dichotomous, but it has to be ensured that they are not defined as factors. \code{D} has to be a
vector, whereas \code{x} and \code{C} have to be matrices. If the model--based approach is used, the colnames of \code{x} and \code{C} have to be used in \code{fmla}!
  Note that the implemented estimates of the 
  attributable risk are only valid if the data has been obtained under a multinomial sampling model!
  }
  
\references{
    Basu S., Landis J.R. (1995)
    Model-based Estimation of Population Attributable Risk under Cross-sectional Sampling.
   \emph{American Journal of Epidemiology}, \bold{142},1338-1343.
   
    Benichou J., Gail M.H. (1989)
    A Delta Method for Implicitly Defined Random Variables.
    \emph{The American Statistician}, \bold{43}, 41-44.
    
    Benichou J. (2001)
    A review of adjusted estimators of attributable risk.
    \emph{Statistical Methods in Medical Research}, \bold{10}, 195-216.
    
    Efron B. (1979)
    Bootstrap methods: another look at the jackknife.
    \emph{Annals of Statistics}, \bold{7},1-26.
    
    Gefeller O. (1992)
    An annotated bibliography on the attributable risk.
    \emph{Biometrical Journal}, \bold{34}, 1007-1012.
    
    Lehnert-Batar A., Pfahlberg A., Gefeller O. (2006)
    Confidence Intervals for Attributable Risk Estimates under Multinomial Sampling.
    \emph{Biometrical Journal}, to appear.
    
    Quenouille M. (1949)
    Approximation tests of correlation in time series.
    \emph{Journal of the Royal Statistical Society, Series B}, \bold{11}, 18-44.
    
    Rubin D.B. (1981)
    The Bayesian Bootstrap.
    \emph{The Annals of Statistics}, \bold{9}, 130-134.
   }
   

\author{Andrea Lehnert-Batar}

\seealso{ \code{\link{PartialAR}}, \code{\link{boot}}}

\examples{
data(icu)
attach(data.frame(icu))

### Computation of crude AR for INF model-free and ###
###  model-based with variance estimates ###
Exp <- matrix(INF,ncol=1)
colnames(Exp) <- "Infection"
AR(D = STA, x = Exp, Var = "delta")
AR(D = STA, x = Exp, model = TRUE, fmla = "STA~Infection", Var="delta")


### Computation of adjusted AR model-free adjusted for SEX and RACE ###
### First coerce variable RACE into dummy matrix as it is not dichotomous! ###
RACENEW <- model.matrix(~as.factor(RACE)-1)
AR(D = STA, x = Exp, C = cbind(SEX,RACENEW), Var="delta", CI="normal") 

### Computation of joint attributable risk for exposure factors INF and TYP ###
### adjusted for RACE and SEX ###
Exp2 <- cbind(INF,TYP)
Conf <- cbind(RACENEW,SEX)
AR(D = STA, x = Exp2, C = Conf, Var="delta")


### Computation of model-based AR for the exposure factor INF, ###\n
### adjusted for TYP, LOC and SEX ###
Conf <- cbind(RACENEW[,-3],SEX)
colnames(Conf) <- c("white","black","sex")
AR(D = STA, x = Exp, C = Conf, model = TRUE, 
    fmla = "STA~Infection+white+black+sex", Var = "delta", CI = "normal")
}
\keyword{ manip }
