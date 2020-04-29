annot_hgcn_FUN <- function(top_t, hgnc=hgnc){
  ## get genes that are reported as they are in the hgnc$symbol column and create their respective tt
  top_t_hgnc <- top_t[(top_t[,1] %in% hgnc$symbol),]
  top_t_hgnc_eIDs <- hgnc$entrez_id[match(top_t_hgnc[,1], hgnc$symbol)]
  top_t_hgnc_df <- data.frame(hGene=top_t_hgnc[,1], hEID=top_t_hgnc_eIDs, top_t_hgnc[,-1], stringsAsFactors = F)

  ## get the remaining genes not reported in the hgnc$symbol column
  top_t_not_hgnc <- top_t[!(top_t[,1] %in% hgnc$symbol),]

  ## get genes that are reported anywhere in HGNC after the hgnc$symbol column filtering steps
  no_cores <- detectCores() - 2
  cl <- makeCluster(no_cores)
  clusterExport(cl,varlist = c("hgnc"))
  hgnc_anywhere <- t(parSapply(cl, top_t_not_hgnc$geneid, FUN = annot_query_FUN,
                               sym = hgnc$symbol, eid = hgnc$entrez_id, syn = hgnc$alias_symbol,
                               psym = hgnc$prev_symbol))
  stopCluster(cl)
  hgnc_anywhere_na <- as.data.frame(hgnc_anywhere) %>% mutate_all(funs(empty_as_na))
  rownames(hgnc_anywhere_na) <- rownames(hgnc_anywhere)
  hgnc_anywhere_na_sel <- hgnc_anywhere_na[!is.na(hgnc_anywhere_na$V1),]
  hgnc_anywhere_df <- data.frame(hGene=hgnc_anywhere_na_sel$V1, hEID=hgnc_anywhere_na_sel$V2,
                                 top_t_not_hgnc[match(rownames(hgnc_anywhere_na_sel),top_t_not_hgnc$geneid),-1],
                                 stringsAsFactors = F)
  ## All remaining genes that are not in HGNC anywhere
  hgnc_nowhere <- top_t_not_hgnc[-(match(rownames(hgnc_anywhere_na_sel),top_t_not_hgnc$geneid)),]
  hgnc_nowhere_df <- data.frame(hGene=hgnc_nowhere$geneid, hEID=NA, hgnc_nowhere[,-1], stringsAsFactors = F)
  ## combined dataframes
  res_all <- do.call("rbind",list(hgnc_anywhere_df, top_t_hgnc_df, hgnc_nowhere_df)) %>% .[order(.$logFC,decreasing = T),]
}
