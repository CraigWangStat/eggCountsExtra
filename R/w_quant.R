# function to compute w with method quant
w_quant <- function(postFEC){

  n <- length(postFEC)
  # find out if there are outliers:
  upperLimit <- quantile(postFEC, probs = 0.75) + 1.5*IQR(postFEC) # outlier definition of boxplot
  truncMean <- mean(postFEC[postFEC < upperLimit])

  quant95upper <- qpois(p=0.95, lambda = truncMean)
  outlier <- (postFEC > quant95upper)

  #weights:
  distPerc <- (postFEC[outlier] - quant95upper)/(max(postFEC) - quant95upper)
  w <- rep(1,n)
  w[outlier] <- 1 - distPerc + 0.01

  wmo <- sum(w[w<1] * postFEC[w<1])/sum(w[w<1])
  postmean <- mean(postFEC)

  return(list(weight = w, wmo = wmo, postmean = postmean))
}

