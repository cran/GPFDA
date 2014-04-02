\name{wrap}
\alias{wrap}

\title{
Wrap the input data into the certain format.
}
\description{
Wrapping function to convert the data format.
}
\usage{
wrap(functional, scale = NULL, testdata = NULL, list = c("traindata", "testdata"), 
  time = NULL, response = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{functional}{
The functional part of the data.
}
  \item{scale}{
The scale part of the data.
}
  \item{testdata}{
Test data.
}
  \item{list}{
Is the wrapping about the training data or testing data? Default to be the first element in the brackets.
}
  \item{time}{
The term gives time points.
}
  \item{response}{
The term gives response.
}
}
\details{
Can wrap the data frame and list into certain format of data.table. For example, see demo('gpfr').
}
\value{
\item{data}{A data table object after wrapping}
}


\author{
Yafeng Cheng
}

\seealso{
  \code{\link{wrap}};  \code{\link{data.table}};  \code{\link{fdapred}};  \code{\link{gpfr}}
}


