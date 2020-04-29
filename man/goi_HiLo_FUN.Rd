\name{goi_HiLo_FUN}
\alias{goi_HiLo_FUN}
\title{
Generate an MAList with over or low expression sample assignment
}
\description{
This functions takes as input a previously saved rds file from the TCGA-GBM Agilent microarrys that we previously processed using the DECL package. For a specific gene of interest, the expression values are median centered and median-z-transformed. Using a predefined upper and lower bound, patients outside these bounds are assigned to the high or low expressing group. Differential gene expression is conducted for the contrast of high-low groups to return a toptable of differentially expressed genes.
}
\usage{
goi_HiLo_FUN(goi, rds_file, l_b = -0.5, u_b = 0.5, rt_fz = 3.2, btb = 10, wd = 6, ht = 6, pix = 600, ...)
}

\arguments{
  \item{goi}{
character, gene of interest symbol}
  \item{rds_file}{
path to the rds file of one of previosuly saved MALists
}
  \item{l_b}{
numeric, lower bound
}
  \item{u_b}{
numeric, upper bound
}
  \item{rt_fz}{
not used
}
  \item{btb}{
not used
}
  \item{wd}{
width argument passed to ggsave
}
  \item{ht}{
height argument passed to ggsave
}
  \item{pix}{
the pixel per inch (ppi) value passed to ggsave
}
  \item{\dots}{
further arguments passed for plotting and file saving parameters of ggplot and ggsave
}
}

\value{
An object of class list containing an MAList of the high/low group assigned patients and the toptable of the high-low contrast.
}

\author{
Ahmed Sadik
}