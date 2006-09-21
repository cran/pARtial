\name{PartialAR}
\alias{PartialAR}
\title{Partial Attributable Risk Estimates, Variances and Confidence Intervals}
\description{
PartialAR derives estimates for partial attributable risks for multiple exposure factors. Variance estimates and confidence intervals
can optionally be returned. Variance estimates are available based on either resampling methods or the delta method. Confidence intervals based on the
BCa method or alternatively constructed with percentiles of the standard normal distribution are available.
}
\usage{
PartialAR(  D, x, w = NULL, model = NULL, fmla, 
            Var = c("none","delta","boot","bayes","jackknife"), 
            CI = c("none","normal","percentile","BCa"), alpha = 0.05, B = 500)
}

\arguments{
  \item{D}{a vector containing a dichotomous indicator variable for the disease status.}
  \item{x}{a matrix containing different exposure variables. Each column of \code{x} indicates the exposure status of one exposure variable,
            which has to be dichotomous. Categorial risk factors have first to be transformed into dummy variables.}
  \item{w}{a weight vector which is
            used to define the resampling technique. If a nonparametric or bayesian bootstrap
            or the jackknife is used, \code{w} can be ignored by the user as it is regulated by
            the input parameter \code{Var}. If else the user wants to use a different
            resampling method, \code{w} can individually be changed.}
  \item{model}{if \code{model=TRUE}, the conditional probabilities necessary for the partial attributable risks are computed by use of estimated
                coefficients from a logistic regression. If \code{model=NULL}, the conditional probabilities are directly estimated from the data.}
  \item{fmla}{if \code{model=TRUE}, \code{fmla} defines the desired form of the logistic regression model and is an obligatory parameter.}
  \item{Var}{a character string indicating the method of variance estimation: \code{Var="delta"} indicates
            a variance estimated derived over the delta method, \code{Var="boot"} and \code{Var="bayes"} means a nonparametric bootstrap and bayesian bootstrap, respectively,
             \code{Var="jackknife"} the Jackknife. \code{Var="none"} is the default!}
  \item{CI}{a character string indicating the method of confidence interval estimation: if \code{CI="normal"} a confidence
            interval constructed by using percentiles from a standard normal distribution is computed. \code{CI="percentile"} and \code{CI="BCa"} yields confidence
            intervals based on the simple percentile and the BCa method, respectively (only possible when \code{Var="boot"} or \code{Var="bayes"}). \code{CI="none"} is the default!}
 \item{alpha}{the probability of error for the estimation of confidence intervals, yielding a $1-\alpha$ confidence level.}
  \item{B}{number of replications for resampling methods.}
}
\details{
The partial attributable risk additively decomposes the joint attributable risk of the risk factors in \code{p} into risk shares for each single
factor. The parameter conceptually corresponds to the Shapley value from cooperative game theory. The partial attributable risk estimates give the possibility
of ranking the risk factors according to their individual relevance for the disease load within the population under study.  
  
If \code{model=NULL}, the conditional probabilities necessary for the computation of the partial attributable risks are directly estimated from the data set.
Else if \code{model=TRUE}, a logistic regression is used for estimating the conditional probabilities. If the model--based approach is used, a formula
defining the form of the logistic regression model has to be given in \code{fmla}. By choosing the model--based approach, the matrix \code{p} must contain
colnames which are used in \code{fmla} (see example!).

All supplied variants of variance estimation yield asymptotic estimates. The variance estimate based on the delta method is computationally least expensive. It is
based on an expansion of the partial attributable risk about its mean by taking a one step Taylor
approximation. If the sample size of the data set is small, the delta method can yield insufficient results. 
The nonparametric bootstrap (\code{Var="boot"}) is based on Efron's bootstrap. Samples are taken out of the underlying data set
and the attributable risk is repeatedly computed for each of those samples. Rubin's bayesian bootstrap (\code{Var="bayes"}) is a weighted method of the bootstrap.
A weight vector is generated from a Dirichlet distribution. For each sampled weight vector, the partial attributable risks are computed with weighted observations.
The jackknife can be seen as an approximation to Efron's nonparametric bootstrap. Replications of the partial attributable risks are computed by leaving out every
single observation once at a time. If the data set contains far more observations as the number of replications normally computed for a bayesian or nonparametric
bootstrap, the computation of a bootstrap is recommended. For the number of replications, \code{B=500} is normally assumed to suffice. The number of replications should
be increased if BCa-intervals are computed. BCa-intervals thus are computationally expensive, but yield very stable results in many situations.

Results of a simulation study showed that the bayesian bootstrap tends to overestimate the variance. In situations where a sufficient amount of observations is available
within the different strata defined by the status of disease and exposures, the use of variance estimates derived over the delta method is advisable, as the computation
is computationally least expensive. If sparse data situations occur, the nonparametric bootstrap combined with the BCa-method 
should be chosen.}
  
\value{
  If only partial attributable risk estimates are computed a vector of length identical to the number of exposure variables in \code{p} is returned.
  If \code{Var=TRUE}, a list containing the vector of partial attributable risks and the vector with corresponding variance estimates is returned.
  If \code{CI=TRUE}, a list containing the vector of partial attributable risks, their variance estimates and a matrix with endpoints of confidence intervals is returned. 
}

\note{The variables in \code{D} and \code{x} have to be dichotomous, but it has to be ensured that they are not defined as factors. \code{D} has to be a
vector, whereas \code{x} has to be a matrix. If the model--based approach is used, the colnames of \code{x} have to be used in \code{fmla}!
  Note that the implemented estimates of the 
  attributable risk are only valid if the data has been obtained under a multinomial sampling model!
  }

\references{
    Cox L A JR. (1985),
    A new measure of attributable risk for public health
    applications.
   \emph{Management Science}, \bold{31},800-813.
    
    Efron B. (1979)
    Bootstrap methods: another look at the jackknife.
    \emph{Annals of Statistics}, \bold{7},1-26.

    Eide G E, Gefeller O. (1995),
    Sequential and average attributable fractions as aids in the selection of preventive strategies.
    \emph{Journal of Clinical Epidemiology}, \bold{48},645-655.
   
    Grömping U, Weimann U. (2004),
    The asymptotic distribution of the partial attributable risk in cross-sectional
    studies.
    \emph{Statistics}, \bold{38},427-438.
    
    Land M, Vogel C, Gefeller O. (2001),
    Partitioning methods for multi-factorial risk attribution.
    \emph{Statistical Methods in Medical Research}, \bold{10}, 217-230.
    
    Quenouille M. (1949)
    Approximation tests of correlation in time series.
    \emph{Journal of the Royal Statistical Society, Series B}, \bold{11}, 18-44.
   }
   

\author{Andrea Lehnert-Batar}

\seealso{ \code{\link{AR}}, \code{\link{boot}}}

\examples{
data(icu)
attach(data.frame(icu))
### Computation of partial attributable risks together with corresponing ###
### confidence intervals based on variance estimates derived from delta method ###
PartialAR(STA,cbind(CAN,INF,TYP,LOC),Var="delta",CI="normal")

### Compuation of partial attributable risks by model-based approach ###
### using a simple main-effects model ###
Exp <- cbind(CAN,INF,TYP,LOC)
colnames(Exp) <- c("cancer","infection","admission","coma")
PartialAR(D = STA, x = Exp, model = TRUE, 
fmla = "STA~cancer+infection+admission+coma", Var = "delta", CI = "normal")

}


\keyword{ manip }
