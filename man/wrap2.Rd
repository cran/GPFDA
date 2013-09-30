\name{wrap2}
\alias{wrap2}

\title{
Wrapping function to get same sample for each repeated curve.
}
\description{
Use over fitted functional object to get the data into matrix with same dimension.
}
\usage{
wrap2(data, new.sample = NULL)
}

\arguments{
  \item{data}{
The object from "wrap" or under certain data.table format.
}
  \item{new.sample}{
The dimension of the new sampled matrix. Default to be the maximum dimension in the data.
}
}
\details{
If the dimension of each repeated curve is same, this step is default to be skipped.
}
\value{
\item{data}{A data.table object}
}
\author{
Yafeng Cheng
}

\seealso{
  \code{\link{wrap}}
}
