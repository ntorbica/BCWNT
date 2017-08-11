# Extract injection sequence from sample names (has to be vector)
getOrd <- function(dat){
  id <- c(rep(0,length(dat)))
  for(i in 1:length(id)){
    id[i] <- as.numeric(gsub('_','',regmatches(dat[i], regexpr('_[0-9]{3}_', dat[i]))))
  }
  return(id)
}
