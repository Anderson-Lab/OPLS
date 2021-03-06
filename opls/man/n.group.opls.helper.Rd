% Generated by roxygen2 (4.0.1.99): do not edit by hand
\name{n.group.opls.helper}
\alias{n.group.opls.helper}
\title{n.group.opls.helper}
\usage{
n.group.opls.helper(X, Y, num_permutations, CV, min_num_OPLS_fact = 0)
}
\arguments{
\item{X}{- n x p matrix, where n is the number of samples and p is the number of variables.}

\item{Y}{- n x 1 matrix. Must be numeric}

\item{num_permutations}{- number of permutation for the randomization test.}

\item{CV}{- number of folds for k-fold cross-validation or -1 for leave one out.}

\item{min_num_OPLS_fact}{- minimum number of OPLS factors to consider.}
}
\description{
An internal helper function. Do not call this directly.
}

