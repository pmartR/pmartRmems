

library(pmartRmems)
library(dplyr)
library(phyloseq)
library(readr)

###############################################################################
# mice data
###############################################################################
e_data <- read.table("inst/extdata/SortedCells_new.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
str(e_data)

f_data <- read.table("inst/extdata/SortedCellsMap.txt", sep = "\t", header = TRUE)
str(f_data)

e_meta <- read.table("inst/extdata/OTUtaxonomy.tsv", sep = "\t", header = TRUE)
str(e_meta)
# colnames for OTUs must match for e_data and e_meta
colnames(e_meta)[1] <- "OTU.ID"

mydata <- as.rRNAdata(e_data, f_data, e_meta,
                      edata_cname = "OTU.ID", fdata_cname = "SampleID", emeta_cname = "OTU.ID")
attributes(mydata)

mice <- mydata
# devtools::use_data(mice)


###############################################################################
# soil data
###############################################################################
e_data <- read.csv("inst/extdata/Incubations_16SFinal.csv", header = TRUE, stringsAsFactors = FALSE)
str(e_data)

load("data-raw/prac_meta_data.RData")
f_data <- prac_meta_data
str(f_data)

# match e_data names with f_data column
names(e_data)
f_data$RNA
names(e_data) <- stringr::str_replace_all(names(e_data), "\\.", "-")

soil <- as.rRNAdata(e_data, f_data, edata_cname = "OTU", fdata_cname = "RNA")
# devtools::use_data(soil)

###############################################################################
# IBD study
###############################################################################

biom <- import_biom("inst/extdata/OTU.biom")

temp <- biom@otu_table@.Data
edata <- temp %>% as.data.frame() %>% mutate(OTU = rownames(temp)) %>% select(OTU, everything())
names(edata) <- stringr::str_replace_all(names(edata), "-", "\\.")

temp <- biom@tax_table@.Data
emeta <- temp %>% as.data.frame() %>% mutate(OTU = rownames(temp)) %>% select(OTU, everything())

temp <- read_delim("inst/extdata/OTU_fdata.txt", "\t",
                   na = c("", "NA", "missing: not provided", "not applicable"))
fdata <- as.data.frame(temp)

ibd <- as.rRNAdata(edata, fdata, emeta, edata_cname = "OTU",
                   fdata_cname = "sample_name", emeta_cname = "OTU")
# devtools::use_data(ibd)
