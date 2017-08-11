# Function to sort the table, primarily by the info
# table, as it seems to contain the injection order
# as part of its description.

align3 <- function(c.ord, to.ord, check.tab, ini){
  # INPUT:
  # c.ord     - The table containing the right order
  # to.ord    - The table to sort by injection order
  # check.tab - The comparison table (input table) to assess for mismatches
  # ini       - initiating index

  # Loop the function call if nrow(to.ord) is exceeded
  # if(ini > nrow(to.ord)){
  #   align2(c.ord, to.ord, check.tab, 1)
  #   # Compare the indexed row, if equal rerun the
  #   # function with next index
  # }

  pcl <- makeCluster(detectCores()-2)
  clusterExport(pcl, c('c.ord', 'to.ord', 'check.tab', 'ini'))

  if(to.ord$SeqNr[ini] == c.ord$SeqNr[ini]){
    if(to.ord[ini,1] == c.ord[ini,1]){
      if(ini == nrow(c.ord)){
        mm <- c(rep(0,nrow(to.ord)))
        cat('Checking for ordering mismatches...\n')
        # for(i in 1:nrow(to.ord)){
        parLapply(pcl, 1:nrow(to.ord), function(i){
          # for(j in 1:nrow(check.tab)){
            lapply(1:nrow(check.tab), function(j){

            if(as.character(to.ord[i,1]) == as.character(check.tab[j,1])){
              # for(k in 2:(ncol(to.ord)-1)){
              lapply(2:(ncol(to.ord)-1), function(k){
                if(to.ord[i,k] != check.tab[j,k]){

                  mm[i] <- mm[i] + 1
                  # to.ord[i,2:ncol(to.ord)] <- check.tab[j,2:ncol(check.tab)]
                }
              })
            }
          })
          gc()
        })
        F.id <- which(mm != 0)

        if(length(F.id != 0)){
          cat('Correcting mismatches...\n')
          # for(i in F.id){
          parLapply(pcl, F.id, function(i){
            to.ord[i,] <- check.tab[which(check.tab[,2] == to.ord[i,2]),]
            gc()
          })
          align3(c.ord, to.ord, check.tab, nrow(c.ord))
        }else{
          stopCluster(pcl)
          return(to.ord)
        }
      }


    }else{
      # print(ini)
      align3(c.ord, to.ord, check.tab, ini+1)
    }

    # If they are inequal, search for the matching
    # row and swap it with the indexed one, then rerun
    # function with the same index
  }else{
    if(ini < (nrow(to.ord)-1)){
      # for(i in (ini+1):nrow(to.ord)){
      parLapply(pcl, (ini+1):nrow(to.ord), function(i){
        if(to.ord$SeqNr[i] == c.ord$SeqNr[ini]){
          if(to.ord[i,1] == c.ord[ini,1]){
            tmp <- to.ord[ini,]
            to.ord[ini,] <- to.ord[i,]
            to.ord[i,] <- tmp
          }
        }
        gc()
      })
    }else{
      # for(i in 1:(nrow(to.ord)-1)){
      parLapply(pcl, 1:(nrow(to.ord)-1), function(i){
        if(to.ord$SeqNr[i] == c.ord$SeqNr[ini]){
          if(to.ord[i,1] == c.ord[ini,1]){
            tmp <- to.ord[ini,]
            to.ord[ini,] <- to.ord[i,]
            to.ord[i,] <- tmp
            # cat(i, 'to', ini,'\n')
          }
        }
        gc()
      })
    }
  }
  align3(c.ord, to.ord, check.tab, ini+1)
}
