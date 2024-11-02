suppressMessages(library(data.table))
suppressMessages(library(optparse))
suppressMessages(library(IRanges))
source("workflow/scripts/version_mapping/s00_mapping_functions.R")

option_list <- list(
  make_option("--input", default=NULL, help="Path and file name of LocusBreaker"),
  make_option("--annot_output", default=NULL, help="Output path and name for cis trans mapping and annotation from Locus Breaker"),
  make_option("--array_path", default=NULL, help="Path to the folder containing the list of targets per array"))
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
LB<-fread(opt$input)
annot_path<-opt$annot_output
array_path<-opt$array_path
####################################
##array version annotation
##load list of targets per array
list_k1<-colnames(fread(paste(array_path,"list_k1.txt",sep="/")))
list_k4<-colnames(fread(paste(array_path,"list_k4.txt",sep="/")))
list_k5<-colnames(fread(paste(array_path,"list_k5.txt",sep="/")))
list_k7<-colnames(fread(paste(array_path,"list_k7.txt",sep="/")))
list_k1_k4_k5_uniprot<-colnames(fread(paste(array_path,"list_k1_k4_k5_uniprot.txt",sep="/")))
##annotate seq id version
annot<-assay_annotation_on_dataset(LB,list_k7,list_k5,list_k4,list_k1)
##annotate uniprot version
annot<-assay_annotation_on_dataset_by_uniprot(annot,list_k1_k4_k5_uniprot)
##save annotated file
fwrite(annot,annot_path,sep=";")
