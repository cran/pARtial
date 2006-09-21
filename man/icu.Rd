\name{icu} 
\alias{icu} 
\non_function{}
\title{ICU Data}
\usage{data(ICU)} \description{ The \code{ICU} data frame
has 200 rows and 7 columns. } \format{
  This data frame contains the following columns:
  \describe{
    \item{STA}{factor, vital status (0 = Lived, 1 = Died).}
    \item{SEX}{factor, sex (0 = male, 1 = female).}
    \item{RACE}{factor, Race (1 = white, 2 = black, 3 = other).}
    \item{CAN}{factor, Cancer part of present problem (0 = No, 1 = Yes).}
    \item{INF}{factor, Infection probable at ICU admission (0 = No, 1 = Yes).}
    \item{TYP}{factor, Type of admission (0 = Elective, 1 = Emergency).}
    \item{LOC}{factor, Level of consciousness at ICU admission (0 = no coma, 1 = deep stupor or coma).}
  }
}

\details{
 The ICU data set consists of a sample of 200 subjects who were part of a much larger study on survival of patients following
 admission to an adult intensive care unit (ICU). The data set presented here contains all variables which significantly explained
 the vital status at hospital discharge. The data were collected at Baystate Medical Center in Springfield, Massachusetts.
 These data are copyrighted by John Wiley \& Sons Inc. and must be acknowledged and used accordingly.
}

\source{
 Hosmer and Lemeshow (2000), Applied Logistic Regression: Second Edition.
 
}

\references{ D.W. Hosmer and S. Lemeshow (2000), \emph{Applied Logistic Regression.} New York: Wiley Series in Probability and Mathematical Statistics }

\keyword{datasets}
