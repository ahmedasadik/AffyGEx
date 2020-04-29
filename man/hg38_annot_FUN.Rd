\name{hg38_annot_FUN}
\alias{hg38_annot_FUN}
\title{
Annotation check using hg38 release
}
\description{
This functions takes as input a toptable and checks both ENST and gene symbols with reference to the hg38 reference assembly. Symbols are updated and entrez gene IDs are added.
}
\usage{
hg38_annot_FUN(top_t, ncbi, hugo, hmart, als_sym)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{top_t}{
toptable
}
  \item{ncbi}{
annotations retreived from the NCBI FTP server
}
  \item{hugo}{
annotations retreived from the HGCN website}
  \item{hmart}{
annotations retreived using biomart}
  \item{als_sym}{
aliases of genes
}
}

\author{
Ahmed Sadik
}