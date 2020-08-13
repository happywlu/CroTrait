# extract O antigen cluster from genome
# OACs 为输入的序列O antigen cluster 序列

O_antigen_clulster <- function(OACs, genome, type = 0) {
  galF <- generate_blast(sequence = OACs, genome = genome, gene = "galF")
  gnd <- generate_blast(sequence = OACs, genome = genome, gene = "gnd")
  
  genome_seq <- tran_name(sequence_file = genome)
  
  if(galF == 0 || gnd == 0) {
    return("bad")
  }
  else{
    galF_split <- unlist(strsplit(galF, "\t"))
    gnd_split <- unlist(strsplit(gnd, "\t"))
    base_range <- c(as.numeric(galF_split[7]),as.numeric(galF_split[8]),as.numeric(gnd_split[7]),as.numeric(gnd_split[8]))
    if(galF_split[1] != gnd_split[1]){
      return("bad")
    }
    else{
      map_galF <- as.numeric(galF_split[4])
      map_gnd <- as.numeric(gnd_split[4])
      len_galF <- contig_length(sequence = OACs, galF_split[2])
      len_gnd <- contig_length(sequence = OACs, gnd_split[2])
      covery_galF <- as.numeric(map_galF)/as.numeric(len_galF)
      covery_gnd <- as.numeric(map_gnd)/as.numeric(len_gnd)
      # 鍒ゆ柇covery鍜宻imilarity
      if(covery_gnd > 0.6 && covery_galF > 0.6 && as.numeric(galF_split[3])> 80 && as.numeric(gnd_split[3]) > 80){
        aim_sequence <- single_seq(sequence = genome_seq, pair_name = galF_split[1])
        # 鍒ゆ柇鏄惁姝ｉ???
        if(as.numeric(galF_split[9]) < as.numeric(galF_split[10]) &&
           as.numeric(gnd_split[9]) < as.numeric(gnd_split[10])) {
          if(type == 0) {
            return("new")
          }
          else if(type == 1) {
            return(as.character(aim_sequence[min(base_range):max(base_range)]))
          }
        }
        # 鍒ゆ柇鏄惁涓鸿礋閾撅紝璐熼摼闇€瑕佸弽鍚戜簰琛?
        else if (as.numeric(galF_split[9]) > as.numeric(galF_split[10]) && 
                 as.numeric(gnd_split[9]) > as.numeric(gnd_split[10])) {
          if(type == 0) {
            return("new")
          }
          else if(type == 1) {
            return(as.character(reverseComplement(aim_sequence[min(base_range):max(base_range)])))
          }
        }
        else {
          return("bad")
        }
      }
      else {
        return("bad")
      }
    }
  }
}

# 
O_serotype <- function(OACs, specie, genome) {
  
  genome_name <- extract_name(genome = genome)
  
  serotypes <- list(c("CO1"), c("DO1a", "DO1b", "DO2"),
                    c("MaO1", "MaO2", "MaO3", "MaO4"), c("MuO1", "MuO2"),
                    c("SO1", "SO2", "SO3", "SO4", "SO6", "SO7"),
                    c("TO1", "TO3", "TO4"), c("UO1"))
  cronobacter <- c("condimenti","dublinensis","malonaticus","muytjensii","sakazakii","turicensis","universalis")
  cronobacter1 <- c("CO", "DO", "MaO", "MuO", "SO", "TO", "UO")
  
  num = 1
  for(i in cronobacter){
    if(specie == i) {
      writeLines(OACs_seq(sequence = OACs, pair_name = c("wzx", "wzm"), specie = cronobacter1[num]), "wzx_wzm.fasta", sep = "\n")
      makeblast_database(infile = "wzx_wzm.fasta", outfile = "wzx_wzm_db")
      blast(query = genome, database = "wzx_wzm_db", outfile = "wzx_wzm.txt")
      if(length(readLines("wzx_wzm.txt")) == 0) {
        new_type <- O_antigen_clulster(OACs = OACs, genome = genome, type = 0)
        return(c(genome_name, specie, new_type))
      }
      else {
        blast_result <- readLines("wzx_wzm.txt", n = 1)
        unlink(c("wzx_wzm_db.nin", "wzx_wzm_db.nhr",
                 "wzx_wzm_db.nsq", "wzx_wzm.fasta","wzx_wzm.txt"), force = TRUE, recursive = FALSE)
        result_split <- unlist(strsplit(blast_result, "\t"))
        len <- contig_length(OACs, result_split[2])
        covery <- as.numeric(result_split[4])/as.numeric(len)
        if(as.numeric(result_split[3]) > 85) {
          if(covery > 0.6) {
            serotype <- unlist(strsplit(result_split[2], "_"))[1]
            return(c(genome_name, specie, serotype))
          }
          else {
            return(c(genome_name, specie, "bad"))
          }
        }
        else {
          new_type1 <- O_antigen_clulster(OACs = OACs, genome = genome, type = 0)
          return(c(genome_name, specie, new_type1))
        }
      }
    }
    num = num + 1
  }
}

# O_serotype(OACs = OACs, specie = "sakazakii", genome = "E.sakazakii.TBCP_3590.fasta")



# MLST <- MLSTs <- tran_name("MLST_seven_genes.fasta")
# ST <- ST_type <- read.table("MLST_type.txt", header = F)
# type == 0 return MLST sequence; type == 1 return ST type
extract_MLST <- function(genome, MLST, ST) {
  seven_genes <- c("atpD","fusA","glnS","gltB","gyrB","infB","pps")
  genome_name <- extract_name(genome = genome)
  mlst <- c()
  for (i in seven_genes) {
    blast_result <- generate_blast(sequence = MLST, genome = genome, gene = i)
    if(blast_result !="bad") {
      blast_split <- unlist(strsplit(blast_result, "\t"))
      if(as.numeric(blast_split[3]) == 100) {
        mlst <- c(mlst, blast_split[2])
      }
      else {return("bad")}
    }
  }
  if(length(mlst) ==7) {
    mlst_seq <- ""
    mlst_type <- ""
    for (n in mlst) {
      n_seq <- single_seq(sequence = MLST, pair_name = n)
      mlst_seq <- paste0(mlst_seq, n_seq)
      mlst_type <- str_c(mlst_type, n, sep = "-")
    }
    
    mlst_type <- str_sub(mlst_type, start = 2)
    mlst_type <- filter(ST, V1 == mlst_type)[1,2]
    
    if (is.na(mlst_type)) {
      return(c(paste0(">", genome_name), mlst_seq, "new"))
    }
    else {
      return(c(paste0(">", genome_name), mlst_seq, mlst_type))
    }
  }
  else {
    return("bad")
  }
}


# construct_tree <- function(genome, genes, MLST = "MLST_seven_genes.fasta", reference = "reference_MLST.fasta",
#                          p1 = "muscle_align_nucleotide.mao", p2 = "infer_ML_nucleotide.mao" )

construct_tree <- function(genome, MLST, ST, reference = "reference_MLST.fasta",
                           p1 = "muscle_align_nucleotide.mao", p2 = "infer_ML_nucleotide.mao" ) {
  
  seven_species <- c("condimenti","dublinensis","malonaticus","muytjensii","sakazakii","turicensis","universalis")
  
  genome_name <- extract_name(genome = genome)
  
  # input reference MLST sequence
  m <- readLines(reference)
  
  # read MLST sequence from genome
  mlst <- extract_MLST(genome = genome, MLST = MLST, ST = ST)
  
  
  if (length(mlst) == 3){
    mlst_seq <- mlst[1:2]
    mlst_type <- mlst[3]
    m <- c(m, mlst_seq)
    writeLines(m, con = "test.fasta", sep = "\n")
    mega(parameter = p1, infile = "test.fasta", outfile = "test")
    mega(parameter = p2, infile = "test.meg", outfile = "test")
    tree <- readLines("test.nwk", n=1)
    unlink(c("test.fasta", "test.meg", "test.nwk", "test_summary.txt"),
           force = TRUE, recursive = FALSE)
    specie <- extract_specie(tree, genome_name)
    if (any(str_detect(seven_species, specie))) {
      return(c(specie, mlst_type))
    }
    else {
      return("bad")
    }
    
  }
  else if (mlst == "bad") {
    return("bad")
  }
}

# construct_tree(genome = "SO2_TBCP-3590.fasta", MLST = MLSTs, ST = ST_type)


# OACs <- tran_name("Cronobacter_OACs.fasta")
# MLSTs <- tran_name("MLST_seven_genes.fasta")
# ST_type <- read.table("MLST_type.txt", header = F)

GermTrait <- function(OACs, genome, ST, MLST, reference = "reference_MLST.fasta",
                      align = "muscle_align_nucleotide.mao", MLtree = "infer_ML_nucleotide.mao"){
  seven_species <- c("condimenti","dublinensis","malonaticus","muytjensii","sakazakii","turicensis","universalis")
  genome_name <- extract_name(genome = genome)
  specie_st <- construct_tree(genome = genome, MLST = MLST, ST = ST, reference = reference, 
                              p1 = align, p2 = MLtree)
  if (length(specie_st) == 2) {
    if(any(str_detect(seven_species, specie_st[1]))) {
      serotype <- O_serotype(OACs = OACs, genome = genome, specie = specie_st[1])
      return(c(serotype, specie_st[2]))
    }
    else {
      return(c(genome_name, "bad", "bad", "bad"))
    }
  }
  else {
    return(c(genome_name, "bad","bad", "bad"))
  }
  
}


GermTrait1 <- function(OACs, genome, ST, MLST, reference = "reference_MLST.fasta",
                       align = "muscle_align_nucleotide.mao", MLtree = "infer_ML_nucleotide.mao") {
  if(file.info(genome)$isdir) {
    my_data <- data.frame(genome = as.character(), species = as.character(), O_type = as.character(), MLST = as.character())
    n=1
    for(i in dir(genome)) {
      m <- GermTrait(OACs = OACs, genome = file.path(genome, i), ST = ST, MLST = MLST, reference = reference,
                     align = align, MLtree = MLtree)
      my_data[n,] <- m
      n = n + 1
      print(n)
    }
    return(my_data)
  }
  else {
    trait <- GermTrait(OACs = OACs, genome = genome, ST = ST, MLST = MLST, reference = reference,
                       align = align, MLtree = MLtree)
    return(trait)
  }
}

