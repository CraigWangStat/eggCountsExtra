# function to compute w with method ratio
w_ratio <- function(preFEC, postFEC){

  n <- length(preFEC)
  ratios <- postFEC/preFEC
  out <- c(which(is.na(ratios) | is.infinite(ratios))) # special case

  #weights:
  w <- rep(1,n)
  w[is.infinite(ratios)] <- 0.01

  for(i in seq(1:n)[-out]){
    if(ratios[i] > 1){ w[i] <- 1/ratios[i] }}

  wmo <- sum(w[w<1] * postFEC[w<1])/sum(w[w<1])
  postmean <- mean(postFEC)

  return(list(weight = w, wmo = wmo, postmean = postmean))
}

