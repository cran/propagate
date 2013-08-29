\name{WelchSatter}
\alias{WelchSatter}
\encoding{latin1}

\title{Welch-Satterthwaite approximation to the 'effective degrees of freedom'}

\description{
Calculates the Welch-Satterthwaite approximation to the 'effective degrees of freedom' by using the samples' uncertainties and degrees of freedoms, as described in Welch (1947) and Satterthwaite (1946).
}

\usage{
WelchSatter(ufinal, usamp, df, alpha = 0.05)
}

\arguments{
  \item{ufinal}{the propagated uncertainty of \eqn{y}.}
  \item{usamp}{the uncertainties of the samples, \eqn{x_i}.}   
  \item{df}{the degrees of freedom of the samples, \eqn{\nu_i}.}  
  \item{alpha}{the significance level for the t-statistic. See 'Details'.}
}

\details{
\deqn{\nu_{\rm{eff}} \approx \frac{u(y)^4}{\sum_{i = 1}^n \frac{u(x_i)^4}{\nu_i}},  \quad k = t(1 - \frac{\alpha}{2}, \nu_{\rm{eff}}), \quad u_{\rm{exp}} = k \cdot u(y)}
}

\value{
A list with the following items:\cr
\item{ws.df}{the 'effective degrees of freedom'.}
\item{k}{the coverage factor for calculating the expanded uncertainty.}
\item{u.exp}{the expanded uncertainty.}
}

\author{
Andrej-Nikolai Spiess
}   

\references{
An Approximate Distribution of Estimates of Variance Components.\cr
Satterthwaite FE.\cr
\emph{Biometrics Bulletin} (1946), \bold{2}: 110-114.

The generalization of "Student's" problem when several different population variances are involved.\cr
Welch BL.\cr
\emph{Biometrika} (1947), \bold{34}: 28-35.
}

\examples{
## Taken from GUM H.1.6, 4).
WelchSatter(32, c(25, 9.7, 2.9, 16.6), c(18, 25.6, 50, 2), alpha = 0.01)
}   

\keyword{matrix}
\keyword{multivariate}
\keyword{algebra}

