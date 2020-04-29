geo_df_fun <- function (x)
{
  geo_br <- readLines(x)
  geo_df <- data.frame()
  Study_names_idx <- grep("^\\d+\\.", geo_br)
  for (i in seq_along(Study_names_idx)) {
    idc <- Study_names_idx[i]
    org_idx <- grep("(Organism:)", geo_br[idc:c(idc + 7)])
    type_idx <- grep("(Type:)", geo_br[idc:c(idc + 7)])
    platform_idx <- grep("(Platform:)", geo_br[idc:c(idc +
                                                       7)], ignore.case = T)
    platforms_idx <- grep("(Platforms:)", geo_br[idc:c(idc +
                                                         7)], ignore.case = T)
    related_idx <- grep("(related platforms)", geo_br[idc:c(idc +
                                                              7)], ignore.case = T)
    accession_idx <- grep("(Accession:)", geo_br[idc:c(idc +
                                                         8)])
    geo_df[i, 1] <- geo_br[idc]
    if (length(geo_br[idc + 1])) {
      geo_df[i, 2] <- geo_br[idc + 1]
    }
    else {
      geo_df[i, 2] <- NA
    }
    if (length(org_idx)) {
      geo_df[i, 3] <- geo_br[idc + org_idx - 1]
    }
    else {
      geo_df[i, 3] <- NA
    }
    if (length(type_idx)) {
      geo_df[i, 4] <- geo_br[idc + type_idx - 1]
    }
    else {
      geo_df[i, 4] <- NA
    }
    if (length(platform_idx)) {
      geo_df[i, 5] <- geo_br[idc + platform_idx - 1]
    }
    else if (length(platforms_idx)) {
      geo_df[i, 5] <- geo_br[idc + platforms_idx - 1]
    }
    else if (length(related_idx)) {
      geo_df[i, 5] <- geo_br[idc + related_idx - 1]
    }
    else {
      geo_df[i, 5] <- NA
    }
    if (length(accession_idx)) {
      geo_df[i, 6] <- geo_br[idc + accession_idx - 1]
    }
    else {
      geo_df[i, 6] <- NA
    }
  }
  num_samples_fun <- function(x) {
    splits <- strsplit(x, split = " ")[[1]]
    num_samples <- splits[length(splits) - 1]
    return(num_samples)
  }
  for (i in seq_along(geo_df[, 5])) {
    if (!is.na(geo_df[i, 5])) {
      geo_df[i, 7] <- num_samples_fun(geo_df[i, 5])
    }
    else {
      geo_df[i, 7] <- NA
    }
  }
  colnames(geo_df) <- c("Study", "Summary", "Organism", "Data type",
                        "Platform", "Accessions", "num of Samples")
  geo_df$Study <- map_chr(geo_df$Study, sub, pattern = "^\\d+\\.",
                          replacement = "")
  geo_df$Summary <- map_chr(geo_df$Summary, sub, pattern = "^\\(Submitter supplied\\)",
                            replacement = "")
  geo_df$Organism <- map_chr(geo_df$Organism, sub, pattern = "(Organism:)",
                             replacement = "")
  geo_df$`Data type` <- map_chr(geo_df$`Data type`, sub, pattern = "(Type:)",
                                replacement = "")
  geo_df$Platform <- map_chr(geo_df$Platform, sub, pattern = "(Platform:)",
                             replacement = "", ignore.case = F)
  geo_df$Platform <- map_chr(geo_df$Platform, sub, pattern = "(Platforms:)",
                             replacement = "", ignore.case = F)
  geo_df$Platform <- map_chr(geo_df$Platform, sub, pattern = "\\s.\\d.*Samples",
                             replacement = "")
  geo_df$Accessions <- map_chr(geo_df$Accessions, gsub, pattern = "(Series)",
                               replacement = "", ignore.case = F) %>% map_chr(., gsub,
                                                                              pattern = "(Accession)", replacement = "", ignore.case = F) %>%
    map_chr(., gsub, pattern = "\t\t", replacement = "",
            ignore.case = F) %>% map_chr(., gsub, pattern = ":",
                                         replacement = "", ignore.case = F)
  geo_df$IDs <- map_chr(geo_df$Accessions, function(x) {
    y <- strsplit(x, "\t")[[1]][2]
  })
  geo_df$IDs <- gsub("ID\\s", "", geo_df$IDs)
  geo_df$Accessions <- map_chr(geo_df$Accessions, function(x) {
    y <- strsplit(x, "\t")[[1]][1]
  })
  geo_df <- geo_df[, c("Study", "num of Samples", "Organism",
                       "Data type", "Platform", "Accessions", "IDs", "Summary")]
  geo_df$`num of Samples` <- as.numeric(geo_df$`num of Samples`)
  geo_df
}
