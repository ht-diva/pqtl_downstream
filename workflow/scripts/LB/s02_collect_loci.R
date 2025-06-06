#!/usr/bin/Rscript


library(data.table)
library(dplyr)

#----------#
# taking variants file as input
args <- commandArgs(trailingOnly = TRUE)
loci_path <- snakemake@input
file_path <- snakemake@output[["ofile"]]
# load param for loci selection
nlp12 <- snakemake@params[["NLP12"]]
mhc <- snakemake@params[["MHC"]]
build <- snakemake@params[["build"]]

# convert input to boolean
nlp12 <- nlp12 %in% c(TRUE, "yes", "true", "TRUE", "Yes", "1")
mhc <- mhc %in% c(TRUE, "yes", "true", "TRUE", "Yes", "1")


#--------------#
# Merge variants file
loci <- tibble(
  rbindlist(
    fill = TRUE,
    lapply(
      loci_path,
      function(x) fread(x, data.table=F, fill = TRUE)
      )
    )
  ) %>%
  arrange(chr) %>%
  filter(!is.na(chr))    #remove trait without significant signals
  #filter(!(chr == 6 & !(end < 28477797 | start > 33448354)))    # remove HLA region

#-------------------------------------#
#        Filter MHC and NLRP12        #
#-------------------------------------#

# define region depending on the genomic build
if (build == "37") {
  # NLP12 gene maps to 54,296,995-54,327,657 in GRCh37, but we use suggested positions by Adam () enlarged by +/-20kb.
  nlp12.start <- 54300000
  nlp12.end   <- 54360000
  hla.start <- 28477797
  hla.end   <- 33448354
  cat("\nTo filter MHC and NLP12 regions, genomic positions are set in build", build, "\n")
} else if (build == "38") {
  # Using liftover.broadinstitute.org resulted in: chr19:53816370-53836078, then expanded it for 20kb
  nlp12.start <- 53796000
  nlp12.end   <- 53856000
  # MHC region maps to chr6:28,510,120-33,480,577 in GRCh38 coordinates.
  hla.start <- 28510120
  hla.end   <- 33480577
  cat("\nTo filter MHC and NLP12 regions, genomic positions are set in build", build, "\n")
}

cat(nrow(loci), "loci-target pair were built")


if (mhc) {
loci <- loci %>%
  arrange(chr) %>%
  #filter(!is.na(chr)) %>%   # remove trait without significant signals
  filter(!(chr == 6 & !(end < hla.start | start > hla.end)))    # remove HLA region
 cat("Removing lead SNPS on MHC region.\n")
}

# remove signals overlapping NLP12 region
if (nlp12) {
  loci <- loci %>% filter(!(chr == 19 & (POS > nlp12.start & POS < nlp12.end)))
  cat("Removing lead SNPs in NLP12 region.\n")
}

if (nlp12|mhc) {
cat(nrow(loci), "remaining after filters.\n")
}

#cat(
  # "\nOf total", nrow(loci),
#  "loci,", nrow(loci) - nrow(ex_mhc),
#  "loci belonging to MHC and", nrow(ex_nlp12) - nrow(ex_mhc),
#  "loci belonging to NLRP12 regions were removed.\n"
#  )


#--------------#
# save the joint results
write.csv(loci, file = file_path, quote = F, row.names = F)
