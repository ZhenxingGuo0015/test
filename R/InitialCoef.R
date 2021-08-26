InitialCoef <- function(x.norm, y.norm,
                        ratio, D, iniby = "lm"){
  #### Initialize Coefficients by different methods
  t1 = Sys.time()
  if(iniby == "lm"){
    R.ini = solve(t(D)%*%D)%*%t(D)%*%mylogit(ratio)
  }else if(iniby == "Beta"){
    #set.seed(12345)
    dat = data.frame(
      y = as.numeric(ratio) +
        runif(length(ratio),
              min=1e-5, max = min(1-ratio)),
      D[, !grepl("Intercept", colnames(D))])
    colnames(dat)[2:ncol(dat)] = colnames(D)[
      !grepl("Intercept", colnames(D))]
    res.fit = betareg(y ~ ., data = dat)
    R.ini = summary(res.fit)$coefficients$mean[, 1]

  }else if(iniby == "Beta-Bin"){
    #### beta-bin regression
    dat = data.frame(y = as.numeric(y.norm),
                     n = as.numeric(x.norm + y.norm),
                     D[, !grepl("Intercept", colnames(D))])
    colnames(dat)[3:ncol(dat)] = colnames(D)[!grepl("Intercept",
                                                    colnames(D))]

    res.fit = betabin(cbind(y, n - y) ~ ., ~ 1, data = dat)
    R.ini = coef(res.fit)
  }
}

