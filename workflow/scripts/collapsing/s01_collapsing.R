suppressMessages(library(data.table))
suppressMessages(library(optparse))
suppressMessages(library(IRanges))
suppressMessages(library(dplyr))

option_list <- list(
  make_option("--input", default=NULL, help="Path and file name of LocusBreaker with lit cis_trans/version/gene and protein name/literature review annotation"),
  make_option("--output", default=NULL, help="Output path and name for collapsed annotated LB"),
  make_option("--mapping", default=NULL, help="Mapping file path for cis and trans"))
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
LB<-fread(opt$input)
output_path<-opt$output
mapping_file<-fread(opt$mapping)


##prepare mapping file
mapping_file$SeqId<-paste("seq.",gsub("-", ".",mapping_file$SeqId),sep="")  
mapping_file$cis_end<-(mapping_file$TSS+500000)
mapping_file$cis_start<-(mapping_file$TSS-500000)


##save original col names and create a locus ID
col<-colnames(LB)
LB$loc_ID<-paste(LB$SNPID,LB$phenotype_id,sep="_")
###divide LB in LB cis and LB trans
LB_cis<-LB  %>%
  filter(cis_or_trans=="cis")
LB_trans<-LB  %>%
    filter(cis_or_trans=="trans")

##fake duplicated file ###TO RM IN SM
#f_LB_trans_2<-rbind(LB_trans[1:10],LB_trans[1:10,])
f_LB_trans_2<-LB_trans



###collapsing trans LB
df_trans<-data.frame(matrix(ncol = ncol(f_LB_trans_2), nrow = 0))
colnames(df_trans) <- colnames(f_LB_trans_2) ###here in SM replace all f_LB_trans_2 by LB_trans
# a<-0
for (i in unique(f_LB_trans_2$loc_ID)){
  # a<-a+1
  # print(a)
  temp<-f_LB_trans_2[f_LB_trans_2$loc_ID==i,]
  out<-temp[1,]
  out$UniProt_ID<-paste(unique(temp$UniProt_ID),collapse="|")
  out$Entrez_Gene_ID<-paste(unique(temp$Entrez_Gene_ID),collapse="|")
  out$Protein.names<-paste(unique(temp$Protein.names),collapse="|")
  out$symbol<-paste(unique(temp$symbol),collapse="|") #check if variables names are the same
  out$UniProt_EntrezID_match<-paste(unique(temp$UniProt_EntrezID_match),collapse="|") 
  out$unip_matching_signal<-paste(unique(temp$unip_matching_signal),collapse="|")
  out$unip_matching_study<-paste(unique(temp$unip_matching_study),collapse="|")
  out$unip_matching_number_ids<-paste(unique(temp$unip_matching_number_ids),collapse="|")
  out$uniprot_match<-ifelse("YES"%in%(temp$uniprot_match),"YES","NO")
  out$new_uniprot<-ifelse("UniProt_previously_assayed"%in%(temp$uniprot_match),"UniProt_previously_assayed","New_UniProt")
df_trans<-rbind(df_trans,out)
}
##restrict to the original colum (in theory not needed - to check)

LB_trans<-df_trans%>%
  select(any_of(col))

##fake duplicated file ###TO RM IN SM
# f_LB_cis_2<-rbind(LB_cis[1:10],f_LB_cis[1:10,])
f_LB_cis_2<-LB_cis
###collapsing cis LB
df_cis<-data.frame(matrix(ncol = ncol(f_LB_cis_2), nrow = 0))
colnames(df_cis) <- colnames(f_LB_cis_2)
a<-0
for (i in unique(f_LB_cis_2$loc_ID)){
   a<-a+1
   print(a)
  temp<-f_LB_cis_2[f_LB_cis_2$loc_ID==i,]
  ir1<-IRanges(start=temp$start[1],end=temp$end[1])
  seqId<-unique(temp$phenotype_id)
  ###define the main row(s)
  map_temp <- mapping_file[mapping_file$target == seqId,]
  for (j in 1:nrow(temp)){
    map_temp_2 <- map_temp[map_temp$symbol == temp$symbol[j], ]
    ir2<-IRanges(start = map_temp_2$cis_start, end = map_temp_2$cis_end)
    ov <- countOverlaps(ir1, ir2)
    temp$MAIN[j]<-ifelse(ov==1,TRUE,FALSE)
  }
  ##create the collapsed row
  out<-temp[1,]
  ##gene name collapse all gene names from main row(s)
  out$Entrez_Gene_ID<-paste(unique(temp$Entrez_Gene_ID[temp$MAIN==TRUE]),collapse="|")
  out$symbol<-paste(unique(temp$symbol[temp$MAIN==TRUE]),collapse="|") #check if variables names are the same
  
  ##uniprot and prot name
  ##check if any unmatch among the main 
  if (TRUE%in%(c(is.na(temp$UniProt_EntrezID_match[temp$MAIN==TRUE]),temp$UniProt_EntrezID_match[temp$MAIN==TRUE]==FALSE))){
    # if yes all uniprot and all prot names 
    out$UniProt_ID<-paste(unique(temp$UniProt_ID),collapse="|")
    out$Protein.names<-paste(unique(temp$Protein.names),collapse="|")
    out$UniProt_EntrezID_match<-paste(unique(temp$UniProt_EntrezID_match),collapse="|") 
    ##version array
    ##check if any unmatch among the main:
    ##if yes, if anny already assessed in all -> "already assessed", otherwise "newly"
    out$new_uniprot<-ifelse("UniProt_previously_assayed"%in%(temp$uniprot_match),"UniProt_previously_assayed","New_UniProt")
    ##lit review uniprot
    ##check if any unmatch among the main:
    ##if yes, if any match it in all -> "match lit", otherwise "new"
    out$unip_matching_signal<-paste(unique(temp$unip_matching_signals),collapse="|")
    out$unip_matching_study<-paste(unique(temp$unip_matching_study),collapse="|")
    out$unip_matching_number_ids<-paste(unique(temp$unip_matching_number_ids),collapse="|")
    out$uniprot_match<-ifelse("YES"%in%(temp$uniprot_match),"YES","NO")
  }else{
    ##if no only collapse the unique uniprot and prot names from main
    out$UniProt_ID<-paste(unique(temp$UniProt_ID[temp$MAIN==TRUE]),collapse="|")
    out$Protein.names<-paste(unique(temp$Protein.names[temp$MAIN==TRUE]),collapse="|")
    out$UniProt_EntrezID_match<-paste(unique(temp$UniProt_EntrezID_match[temp$MAIN==TRUE]),collapse="|") 
     
    ##version array
    ##check if any unmatch among the main:
    ##if no, if any "already assessed" in main-> "already assessed", otherwise "newly"
    out$new_uniprot<-ifelse("UniProt_previously_assayed"%in%(temp$uniprot_match[temp$MAIN==TRUE]),"UniProt_previously_assayed","New_UniProt")
    ##lit review uniprot
    ##check if any unmatch among the main:
    ##if no, if any "match lit" in main-> "match lit", otherwise "new"
    out$unip_matching_signal<-paste(unique(temp$unip_matching_signals),collapse="|")
    out$unip_matching_study<-paste(unique(temp$unip_matching_study),collapse="|")
    out$unip_matching_number_ids<-paste(unique(temp$unip_matching_number_ids),collapse="|")
    out$uniprot_match<-ifelse("YES"%in%(temp$uniprot_match[temp$MAIN==TRUE]),"YES","NO")
      }
 
  df_cis<-rbind(df_cis,out)
}
##restrict to the original colum (in theory not needed - to check)

LB_cis<-df_cis%>%
  select(any_of(col))

LB<-rbind(LB_trans,LB_cis)
LB<-LB[order(chr,POS)]

 
##save
fwrite(LB,output_path)

 
  