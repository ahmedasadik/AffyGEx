\name{TF_rep_final_df_FUN}
\alias{TF_rep_final_df_FUN}
\title{
Define the representation of different transcription factor target genes among differentially regulated genes.
}
\description{
This function saves to disk tables that show the number of target genes of each transcription factor (combined and single motifs) represented in the toptable. Different cutoffs are applied using log2FC, p.value, adj.p.value to show how much change occurs in the number of gene representations after applying these cutoffs. The table includes the log2FC value of the transcription factors. Additionally, a list is returned with both representation tables and transcription factor specific tobtables.
}
\usage{
TF_rep_final_df_FUN(t_table, rep_df, fc_filt, res_path, wd = 6, ht = 6, pix = 300,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{t_table}{
topTable genrated after running the hg38_annot_FUN
}
  \item{rep_df}{
data frame of TF target gene representations
}
  \item{fc_filt}{
numeric, log2 fold change cutoff for selecting transcription factors
}
  \item{res_path}{
path for saving files
}
  \item{wd}{
width argument passed ggsave
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
\author{
Ahmed Sadik}
