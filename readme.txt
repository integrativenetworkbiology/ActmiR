#######################################################################################################
### Code for inferring miRNA activity using IWLS                                                   ####
### written by Eunjee Lee  (eunjee.lee@mssm.edu)                                                   ####
### Input : miRNA ID, miRNA (character). Should be one of IDs of miRNA expression levels           ####
###	    matrix of miRNA expression level (numeric), miRexpr (miRNA by sample)                  ####
###         matrix of mRNA expression level (numeric), expr  (gene by sample)                      ####
###         target gene list for the miRNA (character), target                                     ####
###         cutoff for weighted regression (numeric), cutoff                                       ####
#######################################################################################################

 Usage: miRact <- InfermiRactivity(miRNA, miRexpr, expr, target, cutoff)
   
 Example: 

 source("./Infer_miRactivity_forBioinformatics.R")
 expmat <- read.delim("./data/Expr_dataforBioinformatics.txt", check.names=F, row.names=1)
 miRexpmat <- read.delim("./data/miRexp_dataforBioinformatics.txt", check.names=F, row.names=1)
 targetmat <- read.delim("./data/PredictedTarget_dataforBioinformatics.txt",  check.names=F, row.names=1)
 targets <- rownames(targetmat)[which(!is.na(targetmat[,"MIMAT0000062"]))]
 miRact <- InfermiRactivity("MIMAT0000062", miRexpmat, expmat, targets, 0.35)  
