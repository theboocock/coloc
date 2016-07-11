## Internal function  approx.bf.p.ave
##
## Calculate approximate bayes factors averaging over values.
##
##


approx.bf.p.ave  <- function(p, f, type, N, s, suffix=NULL, average_over_var= c(0.01,.1,0.5){
  if(type=="quant") {
    V <- Var.data(f, N)
  } else {
    V <- Var.data.cc(f, N, s)
  }
  z <- qnorm(0.5 * p, lower.tail = FALSE)
  ## Shrinkage factor: ratio of the prior variance to the total variance
  rs <- average_over / (average_over+ V)
  ## Approximate BF  # I want ln scale to compare in log natural scale with LR diff
  lABF = 0.5 * (log(1-rs) + (rs * z^2))
  lABF = mean(lABF)
  ret <- data.frame(V,z,r,lABF)
  if(!is.null(suffix))
    colnames(ret) <- paste(colnames(ret), suffix, sep=".")
  return(ret)  
}


approx.bf.p.estimates.ave <- function(p,f,type, N,s,suffix=NULL, average_over_var=c(0.01,.01,.5){
  rs <- average_over_var^2/(average_over_var^2 + V)
  lABF = 0.5 * (log(1 - rs) + (rs * z^2))
  lABF = mean(lABF)
  ret <- data.frame(V, z, rs, lABF)
  if(!is.null(suffix))
    colnames(ret) <- paste(colnames(ret), suffix, sep = ".")
  return(ret)
}


region.bf  <-  function(){

}



process.dataset <- function(df1, suffix1,df2,suffix2) {
  message('Processing datasets together because correlation specified')
  nd <- names(d)
  if (! 'type' %in% nd)
    stop("dataset ",suffix,": ",'The variable type must be set, otherwise the Bayes factors cannot be computed')

  if(!(d$type %in% c("quant","cc")))
      stop("dataset ",suffix,": ","type must be quant or cc")
  if(d$type=="cc") {
      if(! "s" %in% nd)
          stop("dataset ",suffix,": ","please give s, proportion of samples who are cases")
      if(! "MAF" %in% nd)
          stop("dataset ",suffix,": ","please give MAF for type cc")
      if(d$s<=0 || d$s>=1)
          stop("dataset ",suffix,": ","s must be between 0 and 1")
  }
  if(d$type=="quant") {
      if(!("MAF" %in% nd || "sdY" %in% nd))
          stop("dataset ",suffix,": ","must give MAF or sdY for type quant")
  }
  
  if("beta" %in% nd && "varbeta" %in% nd && ("MAF" %in% nd || "sdY" %in% nd)) {
    if(length(d$beta) != length(d$varbeta))
      stop("dataset ",suffix,": ","Length of the beta vectors and variance vectors must match")
    if(!("snp" %in% nd))
      d$snp <- sprintf("SNP.%s",1:length(d$beta))
    if(length(d$snp) != length(d$beta))
      stop("dataset ",suffix,": ","Length of snp names and beta vectors must match")
 
    if(d$type == 'quant' & !('sdY' %in% nd)) 
      d$sdY <- sdY.est(d$varbeta, d$MAF, d$N)
    
    df <- approx.bf.estimates(z=d$beta/sqrt(d$varbeta),
                              V=d$varbeta, type=d$type, suffix=suffix, sdY=d$sdY)
    df$snp <- as.character(d$snp)
    return(df)
  }

  if("pvalues" %in% nd & "MAF" %in% nd & "N" %in% nd) {
    if (length(d$pvalues) != length(d$MAF))
      stop('Length of the P-value vectors and MAF vector must match')
    if(d$type=="cc" & !("s" %in% nd))
      stop("Must specify s if type=='cc' and you want to use approximate Bayes Factors")
    if(!("snp" %in% nd))
      d$snp <- sprintf("SNP.%s",1:length(d$pvalues))
    df <- data.frame(pvalues = d$pvalues,
                     MAF = d$MAF,
                     snp=as.character(d$snp))    
    colnames(df)[-3] <- paste(colnames(df)[-3], suffix, sep=".")
    df <- subset(df, df$MAF>0 & df$pvalues>0) # all p values and MAF > 0
    abf <- approx.bf.p(p=df$pvalues, f=df$MAF, type=d$type, N=d$N, s=d$s, suffix=suffix)
    df <- cbind(df, abf)
    return(df)  
  }

  stop("Must give, as a minimum, one of:\n(beta, varbeta, type, sdY)\n(beta, varbeta, type, MAF)\n(pvalues, MAF, N, type)")
}
