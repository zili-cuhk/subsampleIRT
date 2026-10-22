\name{subsampleIRT-package}
\alias{subsampleIRT-package}
\alias{subsampleIRT}
\docType{package}
\title{Target-Specific Optimal Subsampling for Penalized Item Response Models}
\description{
Penalized marginal maximum likelihood for multidimensional, multi-group mixed
2PL/3PL models and target-specific optimal subsampling for large data sets.

A typical analysis:
\preformatted{
dat <- irt_data(Y, group = country, booklet = booklet, weights = w)
mod <- irt_model(dat, dimension = dim, guessing = is_mc)
des <- sub_design(dat, mod, targets = c("difficulty", "joint"))
fit <- sub_fit(des, size = 30000, target = "difficulty")
summary(fit)
}
\code{\link{irt_fit}} fits the full data; \code{\link{sub_size}} chooses the
subsample size for a standard-error tolerance; \code{\link{sub_path}} fits a
sequence of penalties on one subsample.

The evaluation engine for the marginal likelihood is chosen automatically:
dense block-wise products for large group-by-booklet blocks and sparse
products for small ones (typically subsamples). It can be fixed with
\code{options(subsampleIRT.engine = "dense")} or \code{"sparse"}.
}
\author{Zili Liu \email{ZLiu20@mdanderson.org}}
\keyword{package}
