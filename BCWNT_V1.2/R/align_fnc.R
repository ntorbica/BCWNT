# Function to sort the table, primarily by the info
# table, as it seems to contain the injection order
# as part of its description.

align <- function(c.ord, to.ord, ini){
  # INPUT:
  # c.ord     - The table containing the right order
  # to.ord    - The table to sort by injection order
  # ini       - initiating index

  # Compare the indexed row, if equal rerun the
  # function with next index
  if(to.ord$SeqNr[ini] == c.ord$SeqNr[ini]){
    if(ini >= length(c.ord$SeqNr)){
      suppressWarnings(return(to.ord))
    }else{
      ini <- ini+1
      align(c.ord, to.ord, ini)
    }

    # If they are inequal, search for the matching
    # row and swap it with the indexed one, then rerun
    # function with the same index
  }else{
    for(i in ini:nrow(to.ord)){
      if(to.ord$SeqNr[i] == c.ord$SeqNr[ini]){
        if(as.character(to.ord[i,1]) == as.character(c.ord[ini,1])){
          tmp <- to.ord[ini,]
          to.ord[ini,] <- to.ord[i,]
          to.ord[i,] <- tmp
        }
      }
    }
  }
  align(c.ord, to.ord,ini)
}

