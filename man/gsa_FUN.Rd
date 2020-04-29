\name{gsa_FUN}
\alias{gsa_FUN}

\title{
CAMERA Gene Set Analysis
}
\description{
This finction performs gene set test using the method implemented in the camera function of th limma package
}
\usage{
gsa_FUN(top_t = NULL, rma_obj, exprmnt, cntrst = "simple", msig.data.lists, cnt_mat = NULL, d_m = NULL, msigs = NULL, wts = NULL, pltfrm = "AFFY", res_path)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{top_t}{
toptable
}
  \item{rma_obj}{
an rma object or an MAList
}
  \item{exprmnt}{
character vector for the short experiment expression
}
  \item{cntrst}{
character, defaults to simple, which denotes a simple treatment-control design, otherwise pairwise contrasts are used as described in the contrast matrix
}
  \item{msig.data.lists}{
list of MSigDB gene lists
}
  \item{cnt_mat}{
contrast matrix}
  \item{d_m}{
design matrix
}
  \item{msigs}{
only used if platform != affy, numeric for the number of MSigDB gene lists
}
  \item{wts}{
numeric weights of different samples
}
  \item{pltfrm}{
defaults to AFFY, but could be others}
  \item{res_path}{
path to save files
}
}

\author{
Ahmed Sadik
}