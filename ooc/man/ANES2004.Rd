\name{ANES2004}
\alias{ANES2004}
\title{Selected data from the 2004 American National Election Study (ANES)}

\description{This dataframe contains selected survey responses by 1212 respondents to the 2004 American National Election Study (ANES). The first fourteen columns contain issue scale responses and the fifteenth and sixteenth columns list 2004 presidential vote choice and party identification, respectively.}

\usage{data(ANES2004)}

\value{
    \itemize{
    \item{\code{libcon}}{ Seven-point liberal-conservative scale.}
    \item{\code{diplomacy}}{ Seven-point diplomacy-military strength scale.}
    \item{\code{govtspend}}{ Seven-point government spending and services scale.}
    \item{\code{defense}}{ Seven-point defense spending scale.}
    \item{\code{healthinsurance}}{ Seven-point health insurance scale.}
    \item{\code{govtjobs}}{ Seven-point guaranteed jobs and income scale.}
    \item{\code{aidblacks}}{ Seven-point aid to blacks scale.}
    \item{\code{govtfundsabortion}}{ Four-point government funding for abortion scale.}
    \item{\code{partialbirthabortion}}{ Four-point partial-birth abortion scale.}
    \item{\code{environmentjobs}}{ Seven-point environment-jobs scale.}
    \item{\code{deathpenalty}}{ Four-point death penalty scale.}
    \item{\code{gunregulations}}{ Five-point gun regulations scale.}
    \item{\code{womenrole}}{ Seven-point women's role scale.}
    \item{\code{gaymarriage}}{ Three-point gay marriage scale.}
    \item{\code{presvote}}{ 2004 presidential vote choice.}
    \item{\code{partyid}}{ Seven-point party identification scale (strong Democrat to strong Republican).}
    }}

\source{
  American National Election Studies. 2017. \url{http://electionstudies.org}.
}

\references{
Armstrong, David A., Ryan Bakker, Royce Carroll, Christopher Hare, Keith T. Poole, and Howard Rosenthal. 2014. \emph{Analyzing Spatial Models of Choice and Judgment with R}. Boca Raton, FL: CRC Press.
}

\author{Christopher Hare and Keith T. Poole}

\seealso{'\link{ooc}'.}

\examples{
data(ANES2004)
}

\keyword{ideal point estimation, NOMINATE, Bayesian latent variable models}
