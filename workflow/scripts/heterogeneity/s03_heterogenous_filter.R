suppressMessages(library(data.table))
suppressMessages(library(optparse))

option_list <- list(
  make_option("--Isquare_thresh", default=NULL, help="Threshold for I2 heterogeneity statistics"),
  make_option("--input", default=NULL, help="Path and file name of file to filter"),
  make_option("--NEF", default=1, help="Number of effective tests to apply Bonferroni correction"),
  make_option("--het_output", default=NULL, help="Output path and name"))
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
LB<-fread(opt$input)
NEF<-opt$NEF
threshold<-(-log10((5*10^(-8))/as.numeric(NEF)))
I2_threshold<-opt$Isquare_thresh
##the rule will be more generalizable if we do not apply any filtering for significancy
heterogenous<-LB[LB$MLOG10P>=threshold&LB$HETISQ>=I2_threshold,]
fwrite(heterogenous,opt$het_output)
