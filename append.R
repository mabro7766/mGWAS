## set working dir

setwd("/scratch/ngs-bioenergy/25okt2017/")

# wd on windows:
# setwd("X:/groups/biostat/myGWASprojects/RENE/MET")

## load the concatenated file

out = read.delim("output.txt")

head(out)

dim(out)

out$logp= -log10(out$Pval)

## load the identifiers of the traitnrs

ID = read.delim("HTU_average_allID.txt")

head(ID)

dim(ID)

## merge both files
    
merge.out = merge(out, ID, by.x="output_trait1.txt", by.y="outputfile")

write.table(merge.out,file="input_findpeaks.txt",quote=FALSE,col.names=TRUE, row.names=FALSE,sep="\t",na=".")
