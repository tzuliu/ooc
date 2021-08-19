\name{ooc}
\alias{ooc}
\title{Ordered Optimal Classification}

\description{Performs Ordered Optimal Classification, an extension of Poole's (2000) nonparametric unfolding procedure for the analysis of ordinal choice data (e.g., ``Strongly Agree,'' ``Somewhat Agree,'' ``Somewhat Disagree,'' ``Strongly Disagree'') as detailed in Hare, Liu, and Lupton (2018).}

\usage{ooc <- function(votemat, dims=2, minvotes=10, lop=0.0.001, 
	polarity=c(1,1), iter=25, nv.method="svm.reg", cost=1)}

\arguments{
  \item{votemat}{A matrix of ordinal choice data, must be consecutive integers starting with 1. Binary choice should be coded so that 1 = Support, 2 = Oppose. Missing data must be coded as \code{NA}.}
  \item{dims}{Number of dimensions to estimate}
  \item{minvotes}{Minimum number of votes required for a respondent to be included in the analysis}
  \item{lop}{A proportion between 0 and 1, the cut-off used for excluding lopsided votes}
  \item{polarity}{A vector specifying the row number of the respondent(s) constrained to have a positive (i.e., right-wing or conservative) score on each dimension}
  \item{iter}{Number of iterations of the modified Optimal Classification algorithm}
  \item{nv.nethod}{The method used to compute normal vectors at each step of the iteration. Choices include: "oprobit" (Ordered probit regression), "svm.reg" (SVM: regression), "svm.class" (SVM: classification), and "krls" (Kernel regularized least squares)}
  \item{cost}{The cost parameter for SVM (ignored unless "svm.reg" or "svm.class" is used to compute normal vectors).}
}

\value{
A list with the following elements:
\item{respondents}{A matrix containing the respondent estimates}
\item{issues}{A matrix containing the estimates for the (c-1) categories of each issue}
\item{issues.unique}{A matrix containing only the estimates for the unique issues}
\item{fits}{The percentage of all binary choices correctly predicted, the Aggregate Proportional Reduction in Error (APRE) statistic for the binary choices, 
	 the percentage of all ordinal choices correctly predicted, and the Aggregate Proportional Reduction in Error (APRE) statistic for the ordinal choices}
\item{OC.result.dominance.binary}{Standard (binary) OC results for the choice matrix arranged in a binary dominance pattern}
\item{errorcount}{Error counts at each iteration of the Optimal Classification algorithm}
\item{votemat.even.binary}{The matrix of ordinal choices recoded into binary categories that are as balanced as possible}
}

\references{Poole, Keith T. 2000. ``Nonparametric Unfolding of Binary Choice Data.'' \emph{Political Analysis} 8(3): 211--237.

Hare, Christopher, Tzu-Ping Liu, and Robert N. Lupton. 2018. ``What Ordered Optimal Classification Reveals about Ideological Structure, Cleavages, and Polarization in the American Mass Public.'' \emph{Public Choice} 176: 57--78.

Armstrong, David A., Ryan Bakker, Royce Carroll, Christopher Hare, Keith T. Poole, and Howard Rosenthal. 2020. \emph{Analyzing Spatial Models of Choice and Judgment}, second edition. Boca Raton, FL: CRC Press.

Hainmueller, Jens and Hazlett, Chad. 2013. ``Kernel Regularized Least Squares: Reducing Misspecification Bias with a Flexible and Interpretable Machine Learning Approach.'' \emph{Political Analysis} 22(2): 143--168.}

\author{Christopher Hare, Tzu-Ping Liu, and Keith T. Poole}

\examples{
\dontshow{
data(ANES2004)
issuescales <- ANES2004[,1:14]}
\dontrun{
data(ANES2004)
issuescales <- ANES2004[,1:14]
result <- ooc(issuescales, dims=2, minvotes=10, lop=0.001, polarity=c(1,1), iter=25, nv.method="svm.reg", cost=1)
presvote <- as.numeric(ANES2004[,15])
presvote[presvote==2] <- 0
x <- result$respondents[,5]
y <- result$respondents[,6]
plot(x, y, type="n", main="Presidential Vote Choice: 2004")
points(x[presvote==1], y[presvote==1], pch="B", col="red", font=2, cex=1.1)
points(x[presvote==0], y[presvote==0], pch="K", col="blue", font=2, cex=1.1)
}

}

\keyword{ideal point estimation, Optimal Classification, nonparametric statistics}
