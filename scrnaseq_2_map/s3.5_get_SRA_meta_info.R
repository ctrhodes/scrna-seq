source("https://bioconductor.org/biocLite.R")
biocLite("SRAdb")

library(SRAdb)
sqlfile <- getSRAdbFile()
sra_con <- dbConnect(SQLite(),sqlfile)
res <- dbGetQuery(sra_con, "select * from sra_ft where run_accession='SRR1274093'")