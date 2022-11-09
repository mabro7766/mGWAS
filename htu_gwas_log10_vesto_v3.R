#R version 3.1.0 (2014-04-10) -- "Spring Dance"
#Copyright (C) 2014 The R Foundation for Statistical Computing
#Platform: x86_64-w64-mingw32/x64 (64-bit)


### how to use GWAS R script.

t<-as.integer(commandArgs(trailingOnly=TRUE))

cat(t)

### set workingdir: dir where phenotype data is

setwd("/scratch/ngs-bioenergy/")

### load library(qqman)
library(qqman)

### step 1:  you need to source the following scripts
#####################################################

source('/scratch/ngs-bioenergy/emma.r')
## this the emma package from Hyun Min Kang  
## http://mouse.cs.ucla.edu/emma/install.html     

source('/scratch/ngs-bioenergy/amm_gwas_vesto_v3.R')
## this runs the GWAS

# source('/group/biostat/myGWASprojects/ARABIDOPSIS/SCRIPTS/plots_gwas.r')
# use library qqman instead

## step 2: load the genotype and kinship data
##############################################################

load('/scratch/ngs-bioenergy/snps2016.RData')

load('/scratch/ngs-bioenergy/kinship.RData')


## step 3: load your phenotype data and log10 transform
##############################################################

ph <- read.delim("CSPP_TO1_met_t.txt")

logph.tmp = log10(ph[,-1])

ecotype_id = ph[,1]

logph=data.frame(ecotype_id,logph.tmp)

rm(logph.tmp)
rm(ph)

## step 4: run the GWAS analysis
##############################################################

# function call: amm_gwas_vesto_v3<-function(Y,X,K,p=0.001,n=2,run=T,calculate.effect.size=FALSE,
#                                            include.lm=FALSE,use.SNP_INFO=FALSE,report=T) 
# p proportion minor allele freq
# Y phenotype
# X genotype
# K kinship
# n trait column (1st column is ecotype_id)
amm_gwas_vesto_v3(logph,snps2016,kinship,p=0.05,n=t, report=F)

# Sys.time()-a
## the function will output a data.frame with n=number of SNPS rows and 9 columns 

## test data (165 ecotypes, 210K SNPs) run 1.4 mins on vesto's dell latitude , 64-bit version.

outputFilename=paste("output_trait",t-1,".txt",sep="")
write.table(output.sort,file=outputFilename,quote=FALSE,col.names=TRUE, row.names=FALSE,sep="\t",na=".")

gwasResults = output.sort[,c(1,2,3,8)]
colnames(gwasResults) = c("SNP","CHR","BP","P")
jpeg(paste("manhattan_trait",t-1,".jpg",sep=""), width=6, height=5, units="in", res=300, quality=100) 
manhattan(gwasResults)
dev.off()




