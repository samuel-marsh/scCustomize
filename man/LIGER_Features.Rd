% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LIGER_Utilities.R
\name{LIGER_Features}
\alias{LIGER_Features}
\title{Extract Features from LIGER Object}
\usage{
LIGER_Features(liger_object, by_dataset = FALSE)
}
\arguments{
\item{liger_object}{LIGER object name.}

\item{by_dataset}{logical, whether to return list with vector of features for each dataset in
LIGER object or to return single vector of unique features across all datasets in object
(default is FALSE; return vector of unique features)}
}
\value{
vector or list depending on \code{by_dataset} parameter
}
\description{
Extract all unique features from LIGER object
}
\examples{
\dontrun{
# return single vector of all unique features
all_features <- LIGER_Features(liger_object = object, by_dataset = FALSE)

# return list of vectors containing features from each individual dataset in object
dataset_features <- LIGER_Features(liger_object = object, by_dataset = TRUE)
}

}
\concept{liger_object_util}
