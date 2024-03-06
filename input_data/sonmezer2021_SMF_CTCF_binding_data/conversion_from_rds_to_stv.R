
library("GenomicRanges")

binding.sites <- readRDS("~/akitaX1-analyses/input_data/sonmezer2021_SMF_CTCF_binding_data/binding.sites.rds")
binding.sites$rownames <- names(binding.sites)

write.table(x=data.frame(binding.sites), file="./binding.sites.tsv", sep="\t", col.names=TRUE, row.names=TRUE, quote=FALSE)

binding.frequencies <- readRDS("~/akitaX1-analyses/input_data/sonmezer2021_SMF_CTCF_binding_data/binding.frequencies.rds")

write.table(x=data.frame(binding.frequencies), file="./binding.frequencies.tsv", sep="\t", col.names=TRUE, row.names=TRUE, quote=FALSE)
