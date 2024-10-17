map <- function(seqId = character(), chr, start, end, mapping_file) {
  map_temp <- mapping_file[mapping_file$target == seqId, ]
  if (nrow(map_temp) == 0) {
    print("seqId not in mapping file")
    return("trans")
  } else{
    map_temp <- map_temp[map_temp$chromosome == chr, ]
    if (nrow(map_temp) == 0) {
      return("trans")
    } else{
      ir1 <- IRanges(start = map_temp$cis_start, end = map_temp$cis_end)
      ir2 <- IRanges(start = start, end = end)
      ov <- countOverlaps(ir1, ir2)
      if (TRUE %in% ov >= 1) {
        return("cis")
      } else{
        return("trans")
      }
    }
  }
}



assay_annotation_on_dataset <-
  function(lit, list_k7, list_k5, list_k4, list_k1) {
    lit$somamer_version <- rep("new", nrow(lit))
    lit$somamer_version <-
      ifelse(
        lit$phenotype_id %in% list_k5 &
          lit$phenotype_id %in% list_k1,
        "already_in_5k_and_1k",
        lit$somamer_version
      )
    lit$somamer_version <-
      ifelse(
        lit$phenotype_id %in% list_k5 &
          !lit$phenotype_id %in% list_k1,
        "already_in_5k",
        lit$somamer_version
      )
    lit$somamer_version <-
      ifelse(
        lit$phenotype_id %in% list_k1 &
          !lit$phenotype_id %in% list_k5,
        "already_in_1k",
        lit$somamer_version
      )
    lit$somamer_version <-
      ifelse(
        lit$somamer_version == "already_in_5k" &
          lit$phenotype_id %in% list_k4,
        "already_in_5k_and_4k",
        lit$somamer_version
      )
    lit$somamer_version <-
      ifelse(
        lit$somamer_version == "already_in_1k" &
          lit$phenotype_id %in% list_k4,
        "already_in_4k_and_1k",
        lit$somamer_version
      )
    lit$somamer_version <-
      ifelse(
        lit$somamer_version == "already_in_5k_and_1k" &
          lit$phenotype_id %in% list_k4,
        "already_in_5k_4k_and_1k",
        lit$somamer_version
      )
    lit$somamer_version <-
      ifelse(lit$somamer_version== "new" &
               lit$phenotype_id %in% list_k4,
             "already_in_4k",
             lit$somamer_version)
    lit$new_somamer <- ifelse(lit$somamer_version == "new", TRUE, FALSE)
    return(lit)
  }

assay_annotation_on_dataset_by_uniprot <-
  function(lit, list_k1_k4_k5_uniprot) {
    lit$new_uniprot <- ifelse(lit$UniProt_ID$in$list_k1_k4_k5_uniprot,"UniProt_previously_assayed","New_UniProt")
    return(lit)
  }
