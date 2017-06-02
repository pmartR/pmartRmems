
e_data <- read.table("inst/extdata/SortedCells_new.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
str(e_data)

f_data <- read.table("inst/extdata/SortedCellsMap.txt", sep = "\t", header = TRUE)
str(f_data)

e_meta <- read.table("inst/extdata/OTUtaxonomy.tsv", sep = "\t", header = TRUE)
str(e_meta)
# colnames for OTUs must match for e_data and e_meta
colnames(e_meta)[1] <- "OTU.ID"

library(pmartRmems)

mydata <- as.rRNAdata(e_data, f_data, e_meta,
                      edata_cname = "OTU.ID", fdata_cname = "SampleID", emeta_cname = "OTU.ID")
attributes(mydata)


rRNA_obj <- mydata
# devtools::use_data(rRNA_obj)
