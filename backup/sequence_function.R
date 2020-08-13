
# extract single sequence from sequence variable
single_seq <- function(sequence, pair_name){
  for(i in seq(length(sequence))) {
    if(sequence@ranges@NAMES[i] == pair_name) {
      return(sequence[[i]])
    }
  }
}


# extract sequences from sequence variable
extract_seq <- function(sequence, type = 0, pair_name) {
  my_sequence <- c()
  if(type == 0) {
    for(i in seq(length(sequence))) {
        if(str_detect(sequence@ranges@NAMES[i], pair_name)) {
          m1 <- paste0(">", sequence@ranges@NAMES[i])
          m2 <- as.character(single_seq(sequence, sequence@ranges@NAMES[i]))
           my_sequence <- c(my_sequence, m1, m2)
        }
    }
    return(my_sequence)
  }
  
  else if(type == 1) {
    for(i in seq(length(sequence))) {
      name <- sequence@ranges@NAMES[i]
      if(str_detect(name,pair_name) && str_detect(name, "galF")==FALSE && str_detect(name, "gnd")==FALSE) {
        m1 <- paste0(">", sequence@ranges@NAMES[i])
        m2 <- as.character(single_seq(sequence, sequence@ranges@NAMES[i]))
        my_sequence <- c(my_sequence, m1, m2)
      }
    }
    return(my_sequence)
  }
}


# extract certain specie wzx or wzm gene from sequence variable
OACs_seq <- function(sequence, pair_name, specie){
  my_sequence <- c()
  for(i in seq(length(sequence))) {
    name <- sequence@ranges@NAMES[i]
    if(any(str_detect(name,pair_name)) && str_detect(name, specie))  {
      m1 <- paste0(">", sequence@ranges@NAMES[i])
      m2 <- as.character(single_seq(sequence = sequence, sequence@ranges@NAMES[i]))
      my_sequence <- c(my_sequence, m1, m2)
    }
  }
  return(my_sequence)
}



