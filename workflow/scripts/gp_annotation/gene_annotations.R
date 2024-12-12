# Load necessary library
rm(list=ls())

suppressMessages(library(optparse))
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(rtracklayer))
suppressMessages(library(GenomicRanges))

option_list <- list(
  make_option("--input", default=NULL, help="Path and file name of LocusBreaker"),
  make_option("--mapping", default=NULL, help="Mapping file path for cis and trans"),
  make_option("--gtf_file", default=NULL, help="Path to the folder containing the GTF annotation file"),
  make_option("--ncbi_file", default=NULL, help="Path to the folder containing the txt file. It includes the feature table"),
  make_option("--uniprot_file", default=NULL, help="Path to the folder containing the list of uniprotkb IDs"),
  make_option("--output", default=NULL, help="Output path"))
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
lb<-fread(opt$input)
mapping<-fread(opt$mapping)
gtf_file_path<-opt$gtf_file
ncbi<-fread(opt$ncbi_file)
df_uniprot_path<-opt$uniprot_file
output<-opt$output

colnames_merged_with_mapping <- c("chr", "start.x", "end.x", "POS", "SNPID", "EA", "NEA", "EAF", "SEF", "MINF", "MAXF",
                                  "BETA", "SE", "DIRECTION", "MLOG10P", "N", "HETISQ", "HETCHISQ", "HETDF", "LHETP",
                                  "phenotype_id", "cis_or_trans", "UniProt_ID", "SomaScan_UniProt_ID", "Target_Full_Name", 
                                  "Entrez_Gene_ID", "SomaScan_Entrez_Gene_ID", "UniProt_EntrezID_match", "symbol", "Protein.names")

colnames_lb_granges_df <- c("chr", "start.x", "end.x", "POS", "SNPID", "EA", "NEA", "EAF", "SEF", "MINF", "MAXF",
                            "BETA", "SE", "DIRECTION", "MLOG10P", "N", "HETISQ", "HETCHISQ", "HETDF", "LHETP",
                            "phenotype_id", "cis_or_trans", "UniProt_ID", "Protein.names", "SomaScan_UniProt_ID", "Target_Full_Name", 
                            "Entrez_Gene_ID", "SomaScan_Entrez_Gene_ID", "symbol", "UniProt_EntrezID_match", "RSID_merged", "gene", 
                            "description", "gene_tss", "description_tss")

mapping$target<-paste("seq.",gsub("-", ".",mapping$SeqId),sep="")

merged_with_mapping <- lb %>%
  left_join(mapping, by = c("phenotype_id" = "target"), relationship = "many-to-many")

merged_with_mapping <- merged_with_mapping[ ,colnames_merged_with_mapping]

rsids <- fread("/exchange/healthds/pQTL/CHRIS/summary_stats/raw/alias/seq.13530.5.regenie.gz", header = TRUE, sep = "\t")
rsids <- as.data.frame(rsids)
rsids$CHROM <- as.character(rsids$CHROM)

# Join based on chromosome = CHROM and POS = GENPOS
merged_result <- merged_with_mapping %>%
  left_join(rsids %>%
              group_by(CHROM, GENPOS) %>%
              summarize(RSID_merged = paste(RSID, collapse = ","), .groups = 'drop'),
            by = c("chr" = "CHROM", "POS" = "GENPOS"))


lb <- merged_result

gtf <- import(gtf_file_path)
head(gtf)

gtf <- gtf[gtf$type == "gene",]
gtf <- gtf[gtf$gene_biotype == "protein_coding", ]

# Filter for rows where seqid starts with "NC_" (this ensures primary assembly only)
primary_assembly_gtf <- gtf[grep("^NC_", seqnames(gtf)), ]

# Filter out mitochondrial and other unwanted chromosomes (keep only 1-22, X, and Y)
primary_assembly_gtf <- primary_assembly_gtf[grepl("NC_0000(0[1-9]|1[0-9]|2[0-2]|23|24)\\..*", seqnames(primary_assembly_gtf)),]

# Create a new column 'chr' by extracting the chromosome number from the seqid
primary_assembly_gtf$chr <- sub("NC_0000([0-9]{2})\\..*","chr\\1", seqnames(primary_assembly_gtf))

# Convert leading zeros in chromosome numbers to remove them (e.g., chr01 -> chr1)
primary_assembly_gtf$chr <- gsub("chr0","chr", primary_assembly_gtf$chr)

primary_assembly_gtf_tss <- primary_assembly_gtf

# Calculate TSS positions directly while keeping it a GRanges object
tss_positions <- ifelse(strand(primary_assembly_gtf_tss) == "+", start(primary_assembly_gtf_tss), end(primary_assembly_gtf_tss))

# Calculate new start and end ranges
new_start <- pmax(1, tss_positions - 500000)
new_end <- tss_positions + 500000

# Update the ranges in the GRanges object
ranges(primary_assembly_gtf_tss) <- IRanges(start = new_start, end = new_end)


primary_assembly_gtf <- as.data.frame(primary_assembly_gtf)
primary_assembly <- GRanges(
  seqnames = Rle(primary_assembly_gtf$chr),
  ranges = IRanges(start = primary_assembly_gtf$start, end = primary_assembly_gtf$end),
  mcols = primary_assembly_gtf  # This will store all the columns from lb as metadata
)

primary_assembly_gtf_tss <- as.data.frame(primary_assembly_gtf_tss)
primary_assembly_tss <- GRanges(
  seqnames = Rle(primary_assembly_gtf_tss$chr),
  ranges = IRanges(start = primary_assembly_gtf_tss$start, end = primary_assembly_gtf_tss$end),
  mcols = primary_assembly_gtf_tss  # This will store all the columns from lb as metadata
)

lb_granges <- GRanges(
  seqnames = Rle(paste0("chr", lb$chr)),
  ranges = IRanges(start = lb$start, end = lb$end),
  mcols = lb  # This will store all the columns from lb as metadata
)

# Find overlaps between lb_granges and genes_gtf
overlaps <- findOverlaps(lb_granges, primary_assembly)
overlaps_tss <- findOverlaps(lb_granges, primary_assembly_tss)

# Extract the indices of overlapping genes
overlapping_genes <- primary_assembly[subjectHits(overlaps)]
overlapping_genes_tss <- primary_assembly[subjectHits(overlaps_tss)]

# Combine lb_granges with overlapping genes
lb_granges$gene <- NA
lb_granges$description <- NA

lb_granges$gene_tss <- NA
lb_granges$description_tss <- NA

# Assign the overlapping gene IDs to the lb_granges
for (i in seq_along(lb_granges)) {
  overlapping_gene <- unique(overlapping_genes$mcols.gene[queryHits(overlaps) == i])
  lb_granges$gene[i] <- paste(overlapping_gene, collapse = ", ")

  overlapping_gene_tss <- unique(overlapping_genes_tss$mcols.gene[queryHits(overlaps_tss) == i])
  lb_granges$gene_tss[i] <- paste(overlapping_gene_tss, collapse = ", ")

  overlapping_description <- unique(overlapping_genes$mcols.description[queryHits(overlaps) == i])
  lb_granges$description[i] <- paste(overlapping_description, collapse = ", ")

  overlapping_description_tss <- unique(overlapping_genes_tss$mcols.description[queryHits(overlaps_tss) == i])
  lb_granges$description_tss[i] <- paste(overlapping_description_tss, collapse = ", ")
}

# Check the updated lb_granges
head(lb_granges)
lb_granges_df <- as.data.frame(lb_granges)
colnames(lb_granges_df)

colnames(lb_granges_df) <- gsub("^mcols\\.", "", colnames(lb_granges_df))
colnames(lb_granges_df)
lb_granges_df <- lb_granges_df[, colnames_lb_granges_df]
names(lb_granges_df)[names(lb_granges_df) == "start.x"] <- "start"
names(lb_granges_df)[names(lb_granges_df) == "end.x"] <- "end"
names(lb_granges_df)[names(lb_granges_df) == "RSID_merged"] <- "RSID"
names(lb_granges_df)[names(lb_granges_df) == "gene"] <- "genes_overlapping_CDS"
names(lb_granges_df)[names(lb_granges_df) == "description"] <- "descriptions_overlapping_CDS"
names(lb_granges_df)[names(lb_granges_df) == "gene_tss"] <- "genes_overlapping_window_TSS"
names(lb_granges_df)[names(lb_granges_df) == "description_tss"] <- "descriptions_overlapping_window_TSS"

write.table(lb_granges_df, output, sep = "\t", quote = F, row.names = F)
