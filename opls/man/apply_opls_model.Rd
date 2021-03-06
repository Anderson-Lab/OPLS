% Generated by roxygen2 (4.0.1.99): do not edit by hand
\name{apply_opls_model}
\alias{apply_opls_model}
\title{apply_opls_model}
\usage{
apply_opls_model(X, Y, opls_results, new_X)
}
\arguments{
\item{X}{- n x p matrix, where n is the number of samples and p is the number of variables.}

\item{Y}{- n x 1 matrix. Must be numeric}

\item{new_X}{- m x p matrix, where m is the number of samples and p is the number of variables.}

\item{opls_model}{- opls model.}
}
\value{
List containing \item{t}{t-score} \item{t_ortho}{t-orthogonal scores} \item{Y_pred}{the predicted-y values for the samples}
}
\description{
Apply a model to new (or old) data. Computes the t scores and the predicted-y value for each sampel in new_X.
}
\examples{
X <- rand(10,10)
new_X <- rand(5,10)
Y <- rand(10,1)
model <- opls(X,Y,1)
res <- apply_opls_model(X,Y,model,new_X)
}

