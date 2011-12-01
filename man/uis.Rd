\name{uis}
\alias{uis}
\title{UIS Drug Treatment study data}
\description{There are 628 data points in the original data, 575 of which have no missing values.

Variable descriptions:

\tabular{lll}{
Variable  \tab Description                     \tab Codes/Values                   \cr
ID        \tab Identification Code             \tab 1 - 628                        \cr
AGE       \tab Age at Enrollment               \tab Years                          \cr
BECK      \tab Beck DepressionScore            \tab 0.000 - 54.000                 \cr
HC        \tab Heroin/Cocaine Use During       \tab 1 = Heroin & Cocaine           \cr
          \tab 3 Months Prior to Admission     \tab 2 = Heroin Only                \cr 
          \tab                                 \tab 3 = Cocaine Only               \cr
          \tab                                 \tab 4 = Neither Heroin nor Cocaine \cr
IV        \tab History of IV Drug Use          \tab 1 = Never                      \cr 
          \tab                                 \tab 2 = Previous                   \cr
          \tab                                 \tab 3 = Recent                     \cr
NDT       \tab Number of Prior Drug Treatments \tab 0 - 40                         \cr
RACE      \tab Subject's Race                  \tab 0 = White                      \cr
          \tab                                 \tab 1 = Non-White                  \cr
TREAT     \tab Treatment Randomization         \tab 0 = Short                      \cr
          \tab Assignment                      \tab 1 = Long                       \cr
SITE      \tab Treatment Site                  \tab 0 = A                          \cr
          \tab                                 \tab 1 = B                          \cr
LEN.T     \tab Length of Stay in Treatment     \tab Days                           \cr
          \tab (Admission Date to Exit Date)   \tab                                \cr 
TIME      \tab Time to Drug Relapse            \tab Days                           \cr
          \tab (Measured from Admission Date)  \tab                                \cr
CENSOR    \tab Event for Treating Lost to      \tab 1 = Returned to Drugs          \cr
          \tab Follow-Up as Returned to Drugs  \tab     or Lost to Follow-Up       \cr
          \tab                                 \tab 0 = Otherwise                  \cr
Y         \tab log of TIME                     \tab                                \cr  
ND1       \tab Component of NDT                \tab                                \cr  
ND2       \tab Component of NDT                \tab                                \cr
LNDT      \tab                                 \tab                                \cr 
FRAC      \tab Compliance fraction             \tab LEN.T/90 for short trt         \cr
          \tab                                 \tab LEN.T/180 for long trt         \cr
IV3       \tab Recent IV use                   \tab 1 = Yes                        \cr
          \tab                                 \tab 0 = No 

}
}
\usage{data(uis)}
\format{A data frame with dimension 575 by 18.}
\source{Table 1.3 of Hosmer,D.W. and Lemeshow, S. (1998) }
\references{Hosmer,D.W. and Lemeshow, S. (1998) Applied Survival
Analysis: Regression Modeling of Time to Event Data, John Wiley and Sons Inc.,
New York, NY}

\keyword{datasets}
