\name{gsa_bar_plot_FUN}
\alias{gsa_bar_plot_FUN}
\title{
GSA results barplots
}
\description{
This function generates barplots of the top 20 (10 up- and 10 down) regulated enriched gene sets. If the number of up- or down regulated gene sets are less than 10, the remaining gene sets are added from the other group of gene sets. Only gene sets meeting the specific p.value cutoff are plotted, even if the number of gene sets are less than 20.
}
\usage{
gsa_bar_plot_FUN(idx, gsa_ls, ts = 8, res_path, coi, wd, ht, pix = 600)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{idx}{
numeric vector of indeces of gene sets in a gsa result list
}
  \item{gsa_ls}{
gsa result list
}
  \item{ts}{
numeric, text size
}
  \item{res_path}{
path to save files}
  \item{coi}{
cutoff of interest, used for the p.value reported}
  \item{wd}{
width argument passed to ggsave}
  \item{ht}{
height argument passed to ggsave}
  \item{pix}{
the pixel per inch (ppi) value passed to ggsave
}
}
\author{
Ahmed Sadik}