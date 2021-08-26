# conduct inference for DMR calling
DMRInfer <- function(stat, nullModel = "standN"){
  ## calculate pval either by mixture normal, truncated normal or by loc fdr
  if(nullModel == "2mix"){
    ## 2-mixed gaussian with mean1=mean2=0, but sd1 !=sd2
    res.EM = mix.2norm.onlysd(Y = stat[!is.na(stat)], pi = 0.1)
    sd0 = res.EM$sd1
    pval = 2*(1- pnorm(abs(stat), mean = 0, sd = sd0) )
    fdr = p.adjust(pval, method = "fdr")
  }else if(nullModel == "trunN"){
    ## truncated normal
      bounds = seq(0.8, 2, 0.1)
      sd0.range = rep(NA, length(bounds))
      #for (ib in 1:length(bounds)) {
      for (ib in seq_len(length(bounds))) {
        sd0.range[ib] = Uniroot.truncNsd(Y = stat,
                                         a = -bounds[ib],
                                         b = bounds[ib])
        }
      if(max(sd0.range) - min(sd0.range) >= 0.5){
        sd0 = min(sd0.range)
        }else{
          sd0 = sd0.range[length(sd0.range)]
          }
      pval = 2*(1- pnorm(abs(stat), mean = 0, sd = sd0) )
      fdr = p.adjust(pval , method = "fdr")
    }else if(nullModel == "locfdr"){
        #### local fdr
        Loc.fdr = Global.fdr = rep(NA, length(stat))
        idx = (!is.na(stat))
        normstat = stat[idx]
        fdrres = locfdr(zz = normstat, nulltype = 1,plot = 0)
        Loc.fdr[idx] = fdrres$fdr
        #### transfer local fdr to global fdr
        fdr.global = numeric(length(normstat))
        xx = fdrres$mat[, "x"]
        leftbreaks = xx - (xx[2] - xx[1])/2
        leftbreaks[1] = min(normstat)
        rightbreaks = xx + (xx[2] - xx[1])/2
        rightbreaks[length(rightbreaks)] = max(normstat)
        #for (i in 1:length(normstat)) {
        for (i in seq_len(length(normstat))) {
          ind = ((leftbreaks <= (-1) * abs(normstat[i])) |
                   (rightbreaks >= abs(normstat[i])))
          F1l = sum(fdrres$mat[ind, "p1f1"])
          Fl = sum(fdrres$mat[ind, "f"])
          fdr.global[i] = 1 - F1l/Fl
        }
        Global.fdr[idx] =  fdr.global
      }else if(nullModel == "standN"){
        ## standard normal
        pval = 2*(1-pnorm(abs(stat)))
        fdr = p.adjust(pval , method = "fdr")
      }

  if(nullModel == "locfdr") {
    return(res = data.frame(fdr = as.numeric(Global.fdr)))
    }else{
      return(res = data.frame(pvalue = as.numeric(pval),
                              padj = as.numeric(fdr)))
      }
}

