rowCheck <- function(tab, e){
  # pcl <- makeCluster(detectCores()-2)

  s <- 1
  c <- list()
  empty <- c()

  while(s != nrow(tab)){
    check <- tab[s,]
    c[[s]] <- which(check == '')
    if(length(c[[s]]) > e){
      empty <- c(empty, s)
    }
    s <- s + 1
  }
  # parLapply(pcl, 1:nrow(tab), function(s){
  #   check <- tab[s,]
  #   c[[s]] <- which(check == '')
  #   if(length(c[[s]]) > e){
  #     empty <- c(empty, s)
  #   }
  #   gc()
  # })
  # stopCluster(pcl)
  return(empty)
}
