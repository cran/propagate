\name{mixCov}
\alias{mixCov}
\encoding{latin1}

\title{Mixing covariances matrices, raw data, summary data or error values into a single covariance matrix}

\description{This function 'mixes' (aggregates) data from covariances matrices, raw and summary (mu, error) data or single error values into one final covariance matrix suitable for use in \code{\link{propagate}}.
}

\usage{
mixCov(..., use = "everything", method = "pearson")
}

\arguments{
  \item{...}{either covariance matrices, raw data, summary data or error values to be aggregated into a single covariance matrix.}  
  \item{use}{see \code{\link{cov}}.}
  \item{method}{see \code{\link{cov}}.}
}

\details{
'Mixes' (aggregates) data of the following types into a final covariance matrix:\cr
1) matrix/dataframe \eqn{\mathbf{A}} of raw data with \eqn{A_1, A_2, \ldots, A_m} variables and at least two observations per variable. This will be transformed into a covariance matrix.\cr
2) vectors of summary data of type \code{c(mu, sigma)}.\cr
3) covariance matrices \eqn{\mathbf{V}} that are already available.\cr
4) single error values \eqn{\sigma}.\cr

This is accomplished by filling a \eqn{m_1 + m_2 + \ldots + m_n} sized square matrix \eqn{\mathbf{C}} succesively with elements \eqn{1 \ldots m_1, m_1 + 1 \ldots m1 + m2, \ldots, m_n + 1 \ldots m_n + m_{n + 1}} with either covariance matrices at \eqn{C_{m_n + 1 \ldots m_n + m_{n + 1}, m_n + 1 \ldots m_n + m_{n + 1}}} or single variance values on the diagonals at \eqn{C_{m_n, m_n}}.
}

\value{The 'mixed' (aggregated) covariance matrix.}

\author{
Andrej-Nikolai Spiess
}   

\references{
Evaluation of measurement data - Guide to the expression of uncertainty in measurement.\cr
JCGM 100:2008 (GUM 1995 with minor corrections).\cr
\url{http://www.bipm.org/utils/common/documents/jcgm/JCGM_100_2008_E.pdf}.

Evaluation of measurement data - Supplement 1 to the Guide to the expression of uncertainty in measurement - Propagation of distributions using a Monte Carlo Method.\cr
JCGM 101:2008.\cr
\url{http://www.bipm.org/utils/common/documents/jcgm/JCGM_101_2008_E.pdf}.
}

\examples{
#######################################################
## Example in Annex H.4 from the GUM 2008 manual
## (see 'References'), measurement of activity.
## This will give exactly the same values as Table H.8.
data(H.4)
attach(H.4)
T0 <- 60
lambda <- 1.25894E-4
Rx <- ((Cx - Cb)/60) * exp(lambda * tx)
Rs <- ((Cs - Cb)/60) * exp(lambda * ts)

mRx <- mean(Rx)
sRx <- sd(Rx)/sqrt(6)
mRx
sRx

mRs <- mean(Rs)
sRs <- sd(Rs)/sqrt(6)
mRs
sRs

R <- Rx/Rs
mR <- mean(R)
sR <- sd(R)/sqrt(6)
mR
sR

cor(Rx, Rs)

## Definition as in H.4.3.
As <- c(0.1368, 0.0018)
ms <- c(5.0192, 0.005)
mx <- c(5.0571, 0.001)

## We have to scale Rs/Rx by sqrt(6) to get the 
## corresponding covariances.
Rs <- Rs/sqrt(6)
Rx <- Rx/sqrt(6)

## Here we mix the raw data in matrix format 
## and the summary data of the variables.
COV <- mixCov(cbind(Rs, Rx), As, ms, mx)
COV

## Prepare the data for 'propagate'.
MEANS <- c(mRs, mRx, As[1], ms[1], mx[1])
SDS <- c(sRs, sRx, As[2], ms[2], mx[2])
DAT <- rbind(MEANS, SDS)
colnames(DAT) <- c("Rs", "Rx", "As", "ms", "mx")

## This will give exactly the same values as 
## in H.4.3/H.4.3.1.
EXPR <- expression(As * (ms/mx) * (Rx/Rs))
RES <- propagate(EXPR, data = DAT, type = "stat", 
                 use.cov = COV, do.sim = TRUE)

RES$prop
RES$sim
}   

\keyword{array}
\keyword{algebra}
\keyword{multivariate}

