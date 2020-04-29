\name{GO_tts_FUN}
\alias{GO_tts_FUN}
\title{
topGO gene ontology analysis
}
\description{
This function performs gene ontology analysis using the elim algorithm described in the topGO package and a hypergeometric based statistic for enrichment.
}
\usage{
GO_tts_FUN(tt, go_annot, ontlgy, node_size, gn_sel_fun)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tt}{
toptable
}
  \item{go_annot}{
gene ontologies of all genes of the toptable
}
  \item{ontlgy}{
ontology to be tested, can be a value of "BP", "CC" or "MF".
}
  \item{node_size}{
minimum number of genes needed to create a node based on the gene ontologies.
}
  \item{gn_sel_fun}{
function for gene selection
}
}
\author{
Ahmed Sadik
}