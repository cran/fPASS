% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Power_Proj_Test_ufDA.R
\name{Power_Proj_Test_ufDA}
\alias{Power_Proj_Test_ufDA}
\title{Power of the Two-sample Projection-based test for functional data with known (or estimated)
eigencomponents.}
\usage{
Power_Proj_Test_ufDA(
  total_sample_size,
  argvals,
  mean_vector,
  eigen_matrix,
  scores_var1,
  scores_var2,
  weights,
  sig.level = 0.05,
  alloc.ratio = c(1, 1),
  npc_to_pick = ncol(eigen_matrix),
  nsim = 10000
)
}
\arguments{
\item{total_sample_size}{Total sample size combing the two groups, must be a positive integer.}

\item{argvals}{The working grid of timepoints to evaluate the eigenfunctions and the mean functions.
It is preferred to take the working grid as dense grid so that
\eqn{\int [\mu_1(t) - \mu_2(t)]\phi_k(t) \,dt} can be calculated with a required precision.}

\item{mean_vector}{The difference in the mean function evaluated at argvals, must be a numeric vector of length same
as that that of argavls.}

\item{eigen_matrix}{The matrix of eigenfunctions evaluated at argvals,
must be a length(argvals) by K matrix, where K is the number of eigenfunctions.}

\item{scores_var1}{The true (or estimate) of covariance matrix of the shrinkage scores for the first group.
Must be symmetric (\code{is.symmetric(scores_var1) == TRUE}) and positive definite
(\code{chol(scores_var1)} without an error!).}

\item{scores_var2}{The true (or estimate) of covariance matrix of the shrinkage scores for the second group.
Must be symmetric (\code{is.symmetric(scores_var2) == TRUE}) and positive definite
(\code{chol(scores_var2)} without an error!).}

\item{weights}{The weights to put to compute the projection \eqn{\int [\mu_1(t) - \mu_2(t)]\phi_k(t) \,dt},
for each \eqn{k=1,\dots, K}. The integral is numerically approximated as
\code{sum(mean_diff(argvals)*eigen_matrix[,k]*weights)}.}

\item{sig.level}{Significance level of the test, default set at 0.05, must be less than 0.2.}

\item{alloc.ratio}{The allocation ratio of samples in the each group. Note that the eigenfunctions
will still be estimated based on the total sample_size, however, the variance
of the \code{shrinkage} scores (which is required to compute the power function) will be
estimated based on the allocation of the samples in each group. Must be given as vector of
length 2. Default value is set at \code{c(1, 1)}, indicating equal sample size.}

\item{npc_to_pick}{Number of eigenfunction to be used to compute the power. Typically this is
becomes handy when the user want to discard few of the last eigenfunctions,
typically with a very small eigenvalues.}

\item{nsim}{The number of samples to be generated from the alternate distribution of
Hotelling T statistic. Default value is 10000.}
}
\value{
Power of the projection-based test for specified difference in the mean function
and the eigencomponents of the covariance of the functional data.
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}

The function \code{Power_Proj_Test_ufDA()} computes the power of
of the two-sample projection-based test for functional response data
setting, when the group difference, the eigenfunctions of the covariance
of the data are specified at dense grid of time points, along with the
(estimated) covariance of the \code{shrinkage} scores.
}
\details{
The projection-based test first extracts K eigenfunctions from the data, and then
project the mean difference function onto each of the eigenfunctions to obtain a K-dimensional
projection vector that reflects the group difference. Wang (2021) pointed that under the null
hypothesis the covariance of K-dimensional functional principal component analysis (fPCA) scores
are the same, and thus a Hotelling \eqn{T^2} test with assuming equal variance of the shrinkage scores
is a valid test. However, Koner and Luo (2023) pointed out that under the alternate hypothesis,
when the difference is mean is significant, the covariance of the shrinkage scores also differ
between the groups. Therefore, while computing the power of test, we must have to derive the
distribution of the Hotelling \eqn{T^2} statistic under the assumption of unequal variance. The
alogrithm for the power of multivariate Hotelling \eqn{T^2} under unequal variance
is coded in \code{\link[=pHotellingT]{pHotellingT()}} function. This particular function is a wrapper around that
function, which inputs the mean difference as a function, and the eigenfunctions and
the scores, and subsequently call the \code{\link[=pHotellingT]{pHotellingT()}} function to compute the power
under unequal variance. See Koner and Luo (2023) for more details on the
formula of the non-null distribution.
}
\examples{

ngrid          <- 101
interval       <- c(-1,1)
gauss.quad.pts <- gss::gauss.quad(ngrid,interval) # evaluation points
working.grid   <- gauss.quad.pts$pt
mean_fn        <- function(t) {0.4*sin(2*pi*t)}
mean_vector    <- mean_fn(working.grid)
eigen_fn       <- function(t, k){ sqrt(2)*{(k==2)*sin(2*pi*t) + (k==1)*cos(2*pi*t)} }
eigen_matrix   <- cbind(eigen_fn(working.grid,1), eigen_fn(working.grid,2))
mean_proj      <- sapply(1:2, function(r) integrate(function(x)
eigen_fn(x,r)*mean_fn(x), interval[1], interval[2])$value)
sig1           <- diag(2)
sig2           <- 2*diag(2)
alp            <- 0.05
n              <- 100
k              <- ncol(eigen_matrix)
cutoff         <- {(n - 2)*k/(n - k -1)}*qf(1-alp, k, n-k-1)
func_power     <- Power_Proj_Test_ufDA(total_sample_size=n,
argvals=working.grid,
mean_vector = mean_vector, eigen_matrix = eigen_matrix,
scores_var1 = sig1, scores_var2= sig2, weights = gauss.quad.pts$wt,
sig.level=alp, alloc.ratio = c(1,1), npc_to_pick=ncol(eigen_matrix),
nsim = 5e3)

}
\references{
Wang, Qiyao (2021)
\emph{Two-sample inference for sparse functional data,  Electronic Journal of Statistics,
Vol. 15, 1395-1423} \cr
\doi{https://doi.org/10.1214/21-EJS1802}.
}
\seealso{
See \code{\link[=pHotellingT]{pHotellingT()}} and \code{\link[=Sim_HotellingT_unequal_var]{Sim_HotellingT_unequal_var()}} for samples
from Hotelling T distribution.
}
\author{
Salil Koner \cr Maintainer: Salil Koner
\email{salil.koner@duke.edu}
}
