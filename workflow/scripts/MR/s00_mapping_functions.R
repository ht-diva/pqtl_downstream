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
    lit$new_target <- rep("new", nrow(lit))
    lit$new_target <-
      ifelse(
        lit$study_id %in% list_k5 &
          lit$study_id %in% list_k1,
        "already_in_5k_and_1k",
        lit$new_target
      )
    lit$new_target <-
      ifelse(
        lit$study_id %in% list_k5 &
          !lit$study_id %in% list_k1,
        "already_in_5k",
        lit$new_target
      )
    lit$new_target <-
      ifelse(
        lit$study_id %in% list_k1 &
          !lit$study_id %in% list_k5,
        "already_in_1k",
        lit$new_target
      )
    lit$new_target <-
      ifelse(
        lit$new_target == "already_in_5k" &
          lit$study_id %in% list_k4,
        "already_in_5k_and_4k",
        lit$new_target
      )
    lit$new_target <-
      ifelse(
        lit$new_target == "already_in_1k" &
          lit$study_id %in% list_k4,
        "already_in_4k_and_1k",
        lit$new_target
      )
    lit$new_target <-
      ifelse(
        lit$new_target == "already_in_5k_and_1k" &
          lit$study_id %in% list_k4,
        "already_in_5k_4k_and_1k",
        lit$new_target
      )
    lit$new_target <-
      ifelse(lit$new_target == "new" &
               lit$study_id %in% list_k4,
             "already_in_4k",
             lit$new_target)
    lit$new <- ifelse(lit$new_target == "new", TRUE, FALSE)
    return(lit)
  }
