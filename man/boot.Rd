\name{boot}
\alias{boot}
\title{Bootstrap and Jackknife Replications for the Attributable Risk or Partial Attributable Risks}
\description{
  boot computes replications of the (adjusted) attributable risk or partial attributable risks in order to obtain the empirical distribution of the estimated 
  parameter. The computation of replications is either based on the nonparametric bootstrap, the bayesian bootstrap or the jackknife.
}
\usage{
boot(D, x, C = NULL, param = c("AR","PAR"), model = NULL, 
     fmla = fmla, type = c("boot","bayes","jackknife"), B = 500) 
}

\arguments{
  \item{D}{a vector containing a dichotomous indicator variable for the disease status.}
  \item{x}{a matrix containing different exposure variables. Each column of \code{x} indicates the exposure status of one exposure variable,
            which has to be dichotomous. Categorial risk factors have first to be transformed into dummy variables. }
  \item{C}{a matrix containing one or multiple confounding variables for adjusting the attributable risk. Every column of \code{C} must be dichotomous. 
            If \code{C=NULL}, replications of the crude attributable risk for the exposure in \code{x} are computed. This input parameter is only
            valid if \code{param="AR"}.}
  \item{param}{defines, whether replications for the (adjusted) attributable risk (\code{param="AR"}) 
                or for partial attributable risks (\code{param="PAR"}) are computed.}
  \item{model}{if \code{model=TRUE} the model based approach of estimating the (partial) attributable risk is used for the computation of replications.}
  \item{fmla}{if the model--based approach is used, \code{fmla} determines the logistic regression model.}
  \item{B}{Number of replications.}
  \item{type}{A character string indicating the resampling method: \code{type="boot"} computes replications based on nonparametric bootstrap,
                \code{type="jackknife"} uses the jackknife.}

}

\details{The nonparametric bootstrap (\code{type="boot"}) is based on Efron's bootstrap. Samples are taken out of the underlying data set
and the attributable risk is repeatedly computed for each of those samples. Rubin's bayesian bootstrap is a weighted method of the bootstrap.
A weight vector is generated from a Dirichlet distribution. For each sampled weight vector, the attributable risk is computed with weighted observations.
The jackknife can be seen as an approximation to Efron's nonparametric bootstrap. Replications of the attributable risk are computed by leaving out every
singel observation once at a time. If the data set contains far more observations as the number of replications normally computed for a bayesian or nonparametric
bootstrap, the computation of a bootstrap is recommended. For the number of replications, \code{B=500} is usually assumed to suffice. 
}

\value{
  If \code{param="AR"}, a vector of length \code{B} containing replications of the attributable risk is returned.
  Else if \code{param="PAR"}, a matrix with \code{B} rows and \code{ncol(x)} columns with replications of the partial attributable risks is returned.
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
\note{ The variables in \code{D} and \code{x} have to be dichotomous, but it has to be ensured that they are not defined as factors. \code{D} has to be a
vector, whereas \code{x} has to be a matrix. If the model--based approach is used, the colnames of \code{x} have to be used in \code{fmla}!
  Note that the implemented estimates of the 
  attributable risk are only valid if the data has been obtained under a multinomial sampling model!}

 

\seealso{\code{\link{AR}}}

\examples{
data(icu)
attach(data.frame(icu))
### Computation of nonparametric bootstrap replications for the ###
### adjusted AR of INF adjusted for SEX ###
boot(D = STA,x = INF, C = SEX, param = "AR", type="boot")

### Computation of nonparametric bootstrap replications ###
### of partial attributable risks ###
boot(D = STA,x = cbind(CAN,INF,TYP,LOC), param = "PAR", type="boot")
}

\keyword{ manip }
