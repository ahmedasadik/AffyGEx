\name{TF_genes_tt_FUN}
\alias{TF_genes_tt_FUN}
\title{
Return number of TF targets present in tt
}
\description{
This function saves to disk tables that show the number of target genes of each transcription factor (combined and single motifs) represented in the toptable. Different cutoffs are applied using log2FC, p.value, adj.p.value to show how much change occurs in the number of gene representations after applying these cutoffs. The table includes the log2FC value of the transcription factors. Additionally, a list is returned with both representation tables and transcription factor specific tobtables.
}
\usage{
TF_genes_tt_FUN(ptrns, t_table, fc_filt = NA, rtrn_value = "dim", pv = NA, apv = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ptrns}{
vector of gene symbols
}
  \item{t_table}{
topTable genrated after running the hg38_annot_FUN
}
  \item{fc_filt}{
numeric, log2 fold change cutoff
}
  \item{rtrn_value}{
can be "dim" if only the number of genes are required, or "tt_fc" if toptable subsets for TF are required
}
  \item{pv}{
numeric, p value cutoff
}
  \item{apv}{
numeric, adjusted p value cutoff
}
}

\author{
Ahmed Sadik}
