
OACs <- tran_name("Cronobacter_OACs.fasta")

MLSTs <- tran_name("MLST_seven_genes.fasta")

ST_type <- read.table("MLST_typenew.txt", header = F)


# a <- read.table("type_OandMLST.txt", header = T)


# m1 <- filter(a, species == "muytjensii", O_type == "MuO1")
# m2 <- m1$genome
# 
# my_seq <- c()
# for(i in m2) {
#   new_seq <- O_antigen_clulster(OACs = OACs, genome = file.path(getwd(), "genomeFromMLST1", paste0(i, ".fas")), type = 1)
#   my_seq <- c(my_seq, paste0(">", i), new_seq)
# }
# 
# writeLines(my_seq, con = "sakazakii_MuO1.fasta", sep = "\n")


O_antigen_clulster(OACs = OACs, genome = file.path(getwd(), "genomeFromMLST1/id-2714.fas"), type = 1)


O_antigen_clulster(OACs = OACs, genome = "E.sakazakii.TBCP_3590.fasta", type = 1)

GermTrait1(OACs = OACs, genome = "id-25.fas", ST = ST_type, MLST = MLSTs)


wlu <- read.table("type_OandMLST.txt", header = T)



###################

GermTrait(OACs = OACs, genome = "GCA_000017665.1.fasta", ST = ST_type, MLST = MLSTs)

m <- GermTrait1(OACs = OACs, genome = file.path(getwd(), "genome1"), ST = ST_type, MLST = MLSTs)


my_data <- data.frame(genome = as.character(), species = as.character(), O_type = as.character(), ST = as.character())


