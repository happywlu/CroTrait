library("dplyr")
library("stringr")

# execute local-installed makeblastdb program
# example
makeblast_database <- function(makeblastdb = "makeblastdb", infile, outfile, dbtype = "nucl") {
  system2(makeblastdb,
          args = c("-in", infile,
                   "-dbtype", dbtype,
                   "-out", outfile))
}


# execute local-installed blastn program
# example
blast <- function(blast="blastn", query, database, outfmt = 6, outfile = "test.txt") {
  system2(blast,
          args = c("-query", query,
                   "-db", database,
                   "-outfmt", outfmt,
                   "-out", outfile))
}

make_blast <- function(OAC, fnn){
  
  makeblast_database(infile = fnn, outfile = "test_db")
  blast(query = OAC, database = "test_db")
  blast_result <- read.table(file = "test.txt", header = F) %>% 
    distinct(V2, .keep_all = TRUE) %>% 
    mutate(strand = ifelse(V9 < V10, "+", "-")) %>%
    arrange(V7) %>% 
    select(V2, V7, V8, strand)
  unlink(c("test_db.nin", "test_db.nhr", "test_db.nsq", "test.txt"))
  return(blast_result)
}

plot_oac <- function(types, cex = 0.6, srt = 45, color = c("red", "blue", "black")) {
  par(mar=c(0,1,0.5,1))
  
  n=length(types)
  
  if(class(types) != "list") {
    types <- list(types)
  }
  
  my_max <- function(x) {
    de <- c()
    for(i in seq(length(x))) {
      m <- x[[i]]
      de <- c(de, max(m[,3]))
    }
    return(max(de))
  }
  xlim = c(-1000, my_max(types))
  ylim = c(-0.02, 1.02)
  
  plot(x=xlim, y=c(0,0), xlim = xlim, ylim = ylim, type = "n", bty = "n",
       axes = F, xlab = "", ylab = "")
  
  for( i in seq(length(types))){
    if(length(types) <= 8) {
      each_width = 1/8
      inner_width = each_width*0.5
      each_media = (8-i)*each_width + each_width/2
    }
    else {
      each_width = 1/length(types)
      inner_width = each_width*0.5
      each_media = n*each_width - each_width/2
    }
    
    
    oac <- types[[i]]
    type <- strsplit(as.character(oac[1,1]), "_")[[1]][1]
    text(type, x = -800, y = each_media, family = "serif", font = 2, cex = cex*1.2)
    
    
    for(gene in seq(dim(oac)[1])) {
      if(abs(oac[gene,2]-oac[gene,3]) > 400) {
        
        y1 = y2 = each_media - inner_width*0.5*0.6
        y3 = each_media - inner_width*0.5
        y4 = each_media
        y5 = each_media + inner_width*0.5
        y6 = y7 = each_media + inner_width*0.5*0.6
        y = c(y1,y2,y3,y4,y5,y6,y7)
        
        if(oac[gene,4] == "+") {
          x1 = x7 = oac[gene, 2]
          x2 = x3 = x5 = x6 = oac[gene, 3] - 400
          x4 = oac[gene,3]
        }
        else if(oac[gene,4] == "-") {
          x1 = x7 = oac[gene, 3]
          x2 = x3 = x5 = x6 = oac[gene, 2] + 400
          x4 = oac[gene,2]
        }
        x = c(x1,x2,x3,x4,x5,x6,x7)
        
        if(str_detect(oac[gene,1], "galF") || str_detect(oac[gene,1], "gnd")) {
          polygon(x,y,col = color[1], border = NA)
        }
        else if(str_detect(oac[gene,1], "wzx") || str_detect(oac[gene,1], "wzy") ||
                str_detect(oac[gene,1], "wzm") || str_detect(oac[gene,1], "wzt")) {
          polygon(x,y,col = color[2], border = NA)
        }
        else {
          polygon(x,y,col = color[3], border = NA)
        }
        text1 <- strsplit(as.character(oac[gene,1]), "_")[[1]][2]
        text(text1,x = mean(c(oac[gene,2], oac[gene,3])), y=each_media + inner_width/1.1,
             family = "serif", font = 4, cex = cex, srt = srt)
      }
      else if(abs(oac[gene,2]-oac[gene,3]) <= 400) {
        
        y1 = each_media - inner_width*0.5
        y2 = each_media
        y3 = each_media + inner_width*0.5
        y = c(y1,y2,y3)
        
        if(oac[gene,4] == "+") {
          x1 = x3 = oac[gene, 2]
          x2 = oac[gene, 3]
        }
        else if(oac[gene,4] == "-") {
          x1 = x3 = oac[gene, 3]
          x2 = oac[gene, 2]
        }
        x = c(x1,x2,x3)
        
        if(str_detect(oac[gene,1], "galF") || str_detect(oac[gene,1], "gnd")) {
          polygon(x,y,col = color[1], border = NA)
        }
        else if(str_detect(oac[gene,1], "wzx") || str_detect(oac[gene,1], "wzy") ||
                str_detect(oac[gene,1], "wzm") || str_detect(oac[gene,1], "wzt")) {
          polygon(x,y,col = color[2], border = NA)
        }
        else {
          polygon(x,y,col = color[3], border = NA)
        }
        text1 <- strsplit(as.character(oac[gene,1]), "_")[[1]][2]
        text(text1,x = mean(c(oac[gene,2], oac[gene,3])), y=each_media + inner_width/1.1,
             family = "serif", font = 4, cex = cex, srt = srt)
      }
    }
    n = n - 1
  }
}

SO1 <- make_blast(OAC = "SO1.fna", fnn = "SO1.fasta")
SO2 <- make_blast(OAC = "SO2.fna", fnn = "SO2.fasta")
SO3 <- make_blast(OAC = "SO3.fna", fnn = "SO3.fasta")
SO4 <- make_blast(OAC = "SO4.fna", fnn = "SO4.fasta")
SO6 <- make_blast(OAC = "SO6.fna", fnn = "SO6.fasta")
SO7 <- make_blast(OAC = "SO7.fna", fnn = "SO7.fasta")
SO8 <- make_blast(OAC = "SO8.fna", fnn = "SO8.fasta")
plot_oac(list(SO1, SO2, SO3, SO4, SO6, SO7, SO8), cex = 0.65,color = c("#96588a","#d54e21","#464646"))

CO1 <- make_blast(OAC = "CO1.fna", fnn = "CO1.fasta")
plot_oac(CO1, cex = 0.65,color = c("#96588a","#d54e21","#464646"))

MaO1 <- make_blast(OAC = "MaO1.fna", fnn = "MaO1.fasta")
MaO2 <- make_blast(OAC = "MaO2.fna", fnn = "MaO2.fasta")
MaO3 <- make_blast(OAC = "MaO3.fna", fnn = "MaO3.fasta")
MaO4 <- make_blast(OAC = "MaO4.fna", fnn = "MaO4.fasta")
MaO5 <- make_blast(OAC = "MaO5.fna", fnn = "MaO5.fasta")
MaO6 <- make_blast(OAC = "MaO6.fna", fnn = "MaO6.fasta")
plot_oac(list(MaO1,MaO2,MaO3,MaO4,MaO5,MaO6), cex = 0.65,color = c("#96588a","#d54e21","#464646"))



MuO1 <- make_blast(OAC = "MuO1.fna", fnn = "MuO1.fasta")
MuO2 <- make_blast(OAC = "MuO2.fna", fnn = "MuO2.fasta")
MuO3 <- make_blast(OAC = "MuO3.fna", fnn = "MuO3.fasta")
MuO4 <- make_blast(OAC = "MuO4.fna", fnn = "MuO4.fasta")
plot_oac(list(MuO1,MuO2,MuO3,MuO4), cex = 0.65,color = c("#96588a","#d54e21","#464646"))



TO1 <- make_blast(OAC = "TO1.fna", fnn = "TO1.fasta")
TO3 <- make_blast(OAC = "TO3.fna", fnn = "TO3.fasta")
TO4 <- make_blast(OAC = "TO4.fna", fnn = "TO4.fasta")
TO6 <- make_blast(OAC = "TO6.fna", fnn = "TO6.fasta")
TO7 <- make_blast(OAC = "TO7.fna", fnn = "TO7.fasta")
TO8 <- make_blast(OAC = "TO8.fna", fnn = "TO8.fasta")
plot_oac(list(TO1,TO3,TO4,TO6,TO7,TO8), cex = 0.65,color = c("#96588a","#d54e21","#464646"))


UO1 <- make_blast(OAC = "UO1.fna", fnn = "UO1.fasta")
plot_oac(UO1, cex = 0.65,color = c("#96588a","#d54e21","#464646"))


DO1a <- make_blast(OAC = "DO1a.fna", fnn = "DO1a.fasta")
DO1b <- make_blast(OAC = "DO1b.fna", fnn = "DO1b.fasta")
DO2 <- make_blast(OAC = "DO2.fna", fnn = "DO2.fasta")
DO3 <- make_blast(OAC = "DO3.fna", fnn = "DO3.fasta")
DO4 <- make_blast(OAC = "DO4.fna", fnn = "DO4.fasta")
DO5 <- make_blast(OAC = "DO5.fna", fnn = "DO5.fasta")
DO6 <- make_blast(OAC = "DO6.fna", fnn = "DO6.fasta")
plot_oac(list(DO1a,DO1b,DO2,DO3,DO4,DO5,DO6), cex = 0.65,color = c("#96588a","#d54e21","#464646"))

plot_oac(list(CO1,DO1a,DO1b,DO2,DO3,DO4,DO5,DO6,MaO1,MaO2,MaO3,MaO4,MaO5,MaO6,
              MuO1,MuO2,MuO3,MuO4,SO1,SO2,SO3,SO4,SO6,SO7,SO8,
              TO1,TO3,TO4,TO6,TO7,TO8,UO1), cex = 0.5, color = c("#96588a","#d54e21","#464646"))



