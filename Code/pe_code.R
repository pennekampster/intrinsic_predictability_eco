

#' Check time series for some important features
#'
#' @param x A time series of real numbers.
#' @return nothing
#' @examples
#' x <- rnorm(1000)
#' 
Check_time_series <- function(x) {
  n <- length(x)
  print("Warning: Many time series analysis methods assume equally spaced time series, please be careful!")
  if(any(is.na(x)))
    warning("Your time series contains some NAs. Be careful, as some analyses / functions deal well with NAs, and some don't.")
  if(any(duplicated(x)))
    warning("Your time series contains some identical values, if there are many, and if they appear consecutively, this could cause problems.")
  if( sum( (x[1:(n-1)] - x[2:n] )==0) > 0)
    warning("Your time series has some identical consecutive values; this could cause problems.")
  ## check for long series of zeros (same values?)
  NULL
}



#' Check time series for some important features
#'
#' @param x the time series, observed at regular intervals.
#' @param m the number of dimensions to embed x into.
#' @param d the time delay
#' @param as.embed: logical; should we return the embedded time series in the order that embed() would?
#' @return data.frame of embedings
#' 
embedd <- function(x, m, d = 1, indices = FALSE, as.embed = TRUE) {
  n <- length(x) - (m-1)*d
  X <- seq_along(x)
  if(n <= 0)
    stop("Insufficient observations for the requested embedding")
  out <- matrix(rep(X[seq_len(n)], m), ncol = m)
  out[,-1] <- out[,-1, drop = FALSE] + rep(seq_len(m - 1) * d, each = nrow(out))
  if(as.embed)
    out <- out[, seq_len(ncol(out))]
  if(!indices)
    out <- matrix(x[out], ncol = m)
  out
}


# calculate the entropy in bits
entropy <- function(wd){

-sum(wd * log2(wd))

}


#' Get the unweighted distribution of "word"s in a time series.
#'
#' @param x A time series of real numbers.
#' @param word_length The word length, also known as the permutation order.
#' @param tau Time lag (not implemented).
#' @param tie_method The method for dealing with tied values; these are the methods used by the rank()
#' function of the base package. Use "average" to treat ties as the same value,
#' "first" to treat ties as different with first value being treated as larger than second,
#' "last" to treat ties as different with last value being treated as larger than second,
#' "random" to give ties random rank.
#'  @return The word distribution.
#' @examples
#' x <- rnorm(1000)
#' word_distribution(x, 3, 1, order)
# word_distribution <- function(x_emb, tie_method) {
#   words <- apply(x_emb, 1, function(x) paste(rev(rank(x, ties.method=tie_method)), collapse="-"))
#   table(words)/nrow(x_emb)
# }
word_distribution <- function(x_emb, tie_method) {
  words <-  unlist(lapply(lapply(1:nrow(x_emb), function(x) (rev(rank(-x_emb[x,])))), paste, collapse="-"))
  table(words)/nrow(x_emb)
}


#' Weighted distribution of "word"s in a time series.
#'
#' @param x A time series of real numbers.
#' @param word_length The word length, also known as the permutation order.
#' @param tau Time lag (not implemented).
#' @param tie_method The method for dealing with tied values; these are the methods used by the rank()
#' function of the base package. Use "average" to treat ties as the same value,
#' "first" to treat ties as different with first value being treated as larger than second,
#' "last" to treat ties as different with last value being treated as larger than second,
#' "random" to give ties random rank,
#' "noise" to add some noise to the time series to break ties.
#' @param noise_amount How much noise to add to the time series; only used for method "noise".
#' Random numbers from a uniform distribution with maximum 1 * 10^(-noise_amount-1) and minimum zero.
#'  @return The word distribution.
#' @examples
#' x <- rnorm(1000)
#' word_distribution(x, 3, 1, order)
weighted_word_distribution <- function(x_emb, tie_method) {
  words <- apply(x_emb, 1, function(x) paste(rev(rank(x, ties.method = tie_method)), collapse="-"))
  weights <- apply(x_emb, 1, function(x) var(x))
  aggregate(weights, list(words), sum)$x / sum(weights)
  
}


#' Calculate the permuation entropy of a times series
#'
#' @param wd The distribution of counts (e.g., of words).
#' @return The entropy of the distribution.
#' @details Words containing NAs are ignored.
#' @examples
#' x <- rnorm(1000)
#' wd <- word_distribution(x, 3, 1, order)
#' PE(wd)
PE <- function(x, weighted, word_length, tau, tie_method, noise_amount=NA) {
  
  ## check if the time series contains all very very similar values
  if(log10(var(x, na.rm=TRUE)) < -20) {
    warning("Time series has variance < 1e-20, so setting entropy to 0.")
    ent <- 0
  }
  
  if(log10(var(x, na.rm=TRUE)) >= -20) {
    ## Add noise if that method used to prevent ties
    if(tie_method=="noise") {
      tie_method <- "average" ## ties.method used by rank is irrelevant if one adds noise, as there will be no ties
      x <- x + runif(length(x)) * 10^(-noise_amount - 1)
    }
    
    ## Embed the data
    x_emb <- embedd(x=x, m=word_length, d=tau)
    
    ## Remove words containing NAs
    x_emb <- na.omit(x_emb)
    
    
    ## Calculate entropy
    if(!weighted)
      wd <- word_distribution(x_emb, tie_method)
    if(weighted)
      wd <- weighted_word_distribution(x_emb, tie_method)
    
    if(tie_method=="average")
      denom <- log2(2*factorial(word_length))
    if(tie_method!="average")
      denom <- log2(factorial(word_length))
    
    
    ent <- entropy(wd)/denom
  }
  ent
}




#' Entropy of a distribution of counts.
#'
#' @param wd The distribution of counts (e.g., of words).
#' @return The entropy of the distribution.
#' @details Words containing NAs are ignored.
#' @examples
#' x <- rnorm(1000)
#' wd <- word_distribution(x, 3, 1, order)
#' PE(wd)
rolling_PE <- function(x, weighted, word_length, tau, tie_method, noise_amount=NA,
                       windowsize=10){
  
  ## Add noise if that method used to prevent ties
  if(tie_method=="noise") {
    tie_method <- "average" ## ties.method used by rank is irrelevant if one adds noise, as there will be no ties
    x <- x + runif(length(x)) * 10^(-noise_amount - 1)
  }
  
  ## Embed the data
  x_emb <- embedd(x=x, m=word_length, d=tau)
  
  words <- apply(x_emb, 1, function(x) paste(rev(rank(x, ties.method = tie_method, na.last="keep")), collapse="-"))
  weights <- apply(x_emb, 1, function(x) var(x))
  
  
  if(tie_method=="average")
    denom <- log(2*factorial(word_length))
  if(tie_method!="average")
    denom <- log(factorial(word_length))
  
  
  if(!weighted)
    PE <- rollapply(words,
                    width=windowsize,
                    function(x) {
                      ifelse(1-(sum(grepl("NA", x)))/length(x)>threshold_prop_valid_words_in_window,
                             entropy(table(x[!grepl("NA", x)]))/denom,
                             NA)
                    },
                    fill=NA)
  
  
  
  
  if(weighted)
    PE <- rollapply(data.frame(weights, words),
                    width=windowsize,
                    function(z) entropy(aggregate(as.numeric(z[,1]), list(z[,2]), sum)$x / sum(as.numeric(z[,1])))/denom,
                    by.column=FALSE,
                    fill=NA)
  
  
  
  PE
  
}
