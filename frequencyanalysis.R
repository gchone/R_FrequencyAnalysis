
# R functions for frequency analysis
# v0.2.5
# Guénolé Choné - guenole.chone@concordia.ca

# This code is in most part a copy-paste from: 
# http://headwateranalytics.weebly.com/blog/flood-frequency-analysis-in-r
# The code for computing the moment matching estimation of the Pearson type III 
# and log-Pearson type III distributions is taken from the fasstr R library:
# https://github.com/bcgov/fasstr


#' FA_fit Fits a given extreme value distribution to an extreme value series
#' @param series A vector representing an extreme value series (e.g., annual maximum flood)
#' @param distribution A three-character name of the extreme value distribution (see ?dist.list())
#' list of available distribution available at https://www.rdocumentation.org/packages/lmomco/versions/1.4.3/topics/Introduction
#' Log-Pearson type III ('lp3') and log-normal ('lno') added to the list
#' @param method The method used to fit the distribution. Could be "lmom". Can be:
#' "lmom", the L-moments method; "mle", the 
#' maximum likelihood estimation; or "mme", the moment matching estimation for 
#' Pearson type III and log-Pearson type III distributions only.
#' @return A list object with: (1) the input data series, (2) distribution information, 
#' (3, 4, 5, 6, 7) Goodness of fit measurements: chi-square test, AIC, BIC, Kolmogorov-Smirnov test and Anderson-Darling test results
#' @export
#' @import lmomco
#' @import PearsonDS
#' @import ADGofTest
FA_fit <- function(series, distribution, method = "lmom") {
  
  library(lmomco)
  library(PearsonDS)
  library(ADGofTest)
  
  #' Defining an internal function to compute the chi-square test
  #' @param param The distributions, result of lmom2par, mle2par, or fitdistrplus::fitdist
  #' @param data The data used for the FA 
  #' @param mme FALSE if param was produced by lmom2par, mle2par, TRUE if it comes from fitdistrplus::fitdist
  chi_square_test <- function(param, data, mme = FALSE) {
    
    #  number of bins follow Sturges' formula
    bins <- ceiling(log2(length(data)) + 1)
    breaks <- quantile(data, probs = seq(0, 1, length.out = bins + 1))
    
    # Compute observed frequencies
    observed_freq <- hist(data, breaks = breaks, plot = FALSE)$counts
    
    if (!mme) # data comes from lmomco, using L-moments or maximum likelihood
    {
      # lmomco does not have a generic cdf() function, but one for each distribution
      # e.g. cdfpe3(), cdfgev(), etc.
      # Let's dynamiccaly call the right cdf function:
      cdf_function <- get(paste0("cdf", param$type), mode = "function")
      
      # Compute expected probabilities for each bin
      expected_probs <- diff(cdf_function(breaks, param))
      nb_params = length(param$para)
    } else {
      #  LP3 with Moment Matching Estimation
      expected_probs <- diff(PearsonDS::ppearsonIII(breaks, params  = param$estimate))
      nb_params = length(param$estimate)
    }
    
    # Compute expected frequencies
    expected_freq <- expected_probs * length(data)
    
    # Avoid division by zero
    expected_freq[expected_freq == 0] <- 1e-6
    
    # Compute Chi-Square statistic
    chi_sq <- sum((observed_freq - expected_freq)^2 / expected_freq)
    
    # Degrees of freedom: (# bins - # estimated parameters - 1)
    df <- bins - nb_params - 1
    # Compute the p value
    p_value <- 1 - pchisq(chi_sq, df)
    
    return(list(Chi_Square = chi_sq, p_value = p_value, df = df))
  }
  
  distribution <- tolower(distribution)
  transformed <- FALSE
  original_series <- series
  # add log Pearson Type 3 to list of distributions supported
  # by lmomco package
  base.dist <- c('lp3', dist.list())
  # log transform series 
  if( distribution == 'lp3' ) {
    series <- log10(series)
    transformed <- TRUE
    distribution <- 'pe3'
  }
  if( distribution == 'lno' ) {
    series <- log10(series)
    transformed <- TRUE
    distribution <- 'nor'
  }
  
  if (any(distribution %in% base.dist) && (method=="mle" || method=="lmom"))
  {
    if (method=="lmom")
    {
      # compute L-moments
      samLmom <- lmom.ub(series)
      # estimate distribution parameters with the L-moments
      distPar <- lmom2par(samLmom, type = distribution)
      
    }
    else 
    {
      
      # estimate distribution parameters with the maximum likelihood
      init <- lmoms(series)
      init$ratios[3] <- 0 # failure rates for mps and mle are substantially 
      # lowered if starting from the middle of the distribution's 
      # shape to form the initial parameters for init.para
      distPar <- mle2par(series, type=distribution, init.lmr=init) # method of MLE

    }
  
    
    if (!is.null(distPar))
    {
      
      # Compute AIC 
      if (transformed) 
      {
        # loglike needs to account for the transformation
        loglike = sum(log(dlmomco(series, distPar)) - log(10^series*log(10)))
      }
      else
      {
        loglike = sum(log(dlmomco(series, distPar)))
      }
      AIC <- 2*length(distPar$para) - 2*loglike
      BIC <- log(length(series))*length(distPar$para) - 2*loglike
      

      # Kolmogorov-Smirnov tests
      cdf_function <- get(paste0("cdf", distPar$type), mode = "function") # find the right cdf function
      ks_result <- ks.test(series, function(x) cdf_function(x, distPar))
      ad_result <- ADGofTest::ad.test(series, cdf_function, distPar)
      
      chi_square <- chi_square_test(distPar, series)
      
    
      results <- list(
        data = original_series,
        method = method,
        distribution = list(
          name = distribution,
          logTransformed = transformed,
          parameters = distPar),
        chi_square = chi_square,
        AIC = AIC,
        BIC = BIC,
        KS = ks_result,
        AD  = ad_result
        ) 
    }
    
  }
  else if ((distribution == "lp3" || distribution == 'pe3' ) && method=="mme")
  {
    # Adding the Moment Matching Estimation for the Pearson type III
    if( distribution == 'pe3' || distribution == 'lp3') {
      ePIII <- function(x, order){
        # compute (centered) empirical centered moments of the data
        if(order == 1) return(mean(x))
        if(order == 2) return(stats::var(x))
        if(order == 3) return(e1071::skewness(x, type = 2))
      }
      dPIII <<- function(x, shape, location, scale) PearsonDS::dpearsonIII(x, shape, location, scale, log = FALSE)
      pPIII <<- function(q, shape, location, scale) PearsonDS::ppearsonIII(q, shape, location, scale, lower.tail = TRUE, log.p = FALSE)
      qPIII <<- function(p, shape, location, scale) PearsonDS::qpearsonIII(p, shape, location, scale, lower.tail = TRUE, log.p = FALSE)
      mPIII <<- function(order, shape, location, scale){
        # compute the empirical first 3 moments of the PIII distribution
        if(order == 1) return(location + shape * scale)
        if(order == 2) return(scale * scale * shape)
        if(order == 3) return(2 / sqrt(shape) * sign(scale))
      }
      # Note that the above forgot to mulitply the scale by the sign of skewness .
      # Refer to Page 24 of the Bulletin 17c
      m <- mean(series)
      v <- stats::var(series)
      s <- stats::sd(series)
      g <- e1071::skewness(series, type = 2)
      
      # This can be corrected, but HEC Bulletin 17b does not do these corrections
      # Correct the sample skew for bias using the recommendation of
      # Bobee, B. and R. Robitaille (1977). "The use of the Pearson Type 3 and Log Pearson Type 3 distributions revisited."
      # Water Resources Reseach 13(2): 427-443, as used by Kite
      #n <- length(series)
      #g <- g*(sqrt(n*(n-1))/(n-2))*(1+8.5/n)
      # We will use method of moment estimates as starting Values for the MLE search
      my_shape <- (2 / g) ^ 2
      my_scale <- sqrt(v) / sqrt(my_shape) * sign(g)
      my_location <- m - my_scale * my_shape
      
      
      start <- list(shape = my_shape, location = my_location, scale = my_scale)
      fit_mme <- fitdistrplus::fitdist(series, 'PIII', start = start, calcvcov = FALSE,
                                   method = "mme", order = 1:3, memp = ePIII, control = list(maxit = 1000))
      
      
      chi_square <- chi_square_test(fit_mme, series, mme = TRUE)
      
      if( distribution == 'pe3' || distribution == 'lp3') {
        results <- list( 
          data = original_series,
          method = method,
          distribution = list(
            name = distribution,
            logTransformed = transformed,
            parameters = fit_mme),
          chi_square = chi_square
        )
      }
      
    }
    return(results)
  } else {
    stop(
      sprintf('Distribution or method not recognized!'))
  }

}


#' FA_getvalues Get flow values for given return periods 
#' @param fa The result of FA_fit
#' @param nep A vector of non-exceedance probabilities. Default values are for the 2, 20, 100 and 350-year return period
#' @param ci The confidence interval to compute. 95% by default. Not used if BootstrapCI was not run on the FA 
#' @return A dataframe the non-exceedance probabilities, their corresponding return period, and the computed flow.
#' if fa is the result of BootstrapCI, the returned dataframe also include the flow for the confidence interval
FA_getvalues <- function(fa, nep = c(0.5, 0.95, 0.99, 1-1/350), ci=0.95) {
  if (fa$method == "lmom" || fa$method == "mle")
  {
    # compute quantiles 
    quant <- par2qua(f = nep, para = fa$distribution$parameters)
  }
  else if (fa$method == "mme")
  {
    quant <- stats::quantile(fa$distribution$parameters, prob = nep)
    quant <- unlist(quant$quantiles)
  }
  if(fa$distribution$logTransformed) {
    quant <- 10^quant
  }
  output = data.frame(nep = nep, rp = prob2T(nep), estimate = quant)
  
  if (!is.null(fa$ci.parameters))
  {
    n.resamples <- dim(fa$ci.parameters)[1]
    quantile.estimates <- matrix(NA, nrow = n.resamples, ncol = length(nep), 
                                 dimnames = list(NULL, nep) ) 
    if (!fa$ci.bootstrap == "mme")
    {
      # begin bootstrapping procedure
      for(i in 1:n.resamples) 
      {

        para <- list(para=fa$ci.parameters[i,], type=fa$distribution$name)
        # estimate quantiles at NEP
        estimated <- qlmomco(nep, para)
        
        # convert quantile estimates to real values if
        # distribution was transformed
        if(fa$distribution$logTransformed) 
        {
          estimated <- 10^estimated
          
        }
        
        # store the quantiles at the desired AEP values
        quantile.estimates[i,] <- estimated
      }

    } else
    {
      fit_mme <- fa$distribution$parameters
      # begin bootstrapping procedure
      for(i in 1:n.resamples) 
      {
        fit_mme$estimate <- fa$ci.parameters[i,]
        
        quant_mme <- stats::quantile(fit_mme, prob = nep)
        quant_mme <- unlist(quant_mme$quantiles)
        if(fa$distribution$logTransformed)
          {quant_mme <- 10 ^ quant_mme}
        
        # store the quantiles at the desired AEP values
        quantile.estimates[i,] <- quant_mme
      }
    }
    # now calculate confidence limits for quantiles
    p <- c((1-ci)/2, (1+ci)/2)
    ci <- sapply(colnames(quantile.estimates), 
                 FUN=function(x){
                   quantile(quantile.estimates[,x], probs=p, na.rm=TRUE)})

    output = data.frame(nep = nep, rp = prob2T(nep),
                        lower = as.vector(ci[1,]), 
                        estimate = quant,
                        upper = as.vector(ci[2,]))
    
  }
  
  return (output)
}




#' BootstrapCI Conducts bootstrap to randomly sample an extreme value series 'n' times for a 
#'  specified distribution to estimate confidence interval for each given 
#'  non-exceedance probability.
#' @param fa Fitted distribution from FA_fit
#' @param n.resamples An integer representing number of re-samples to conduct
#' @param use_lmom TRUE by default, meaning that the confidence intervals are computed using the L-moments method,
#' even when the fitted distribution fa was computed with the maximum likelihood estimation or the moment matching estimation
#' if FALSE, use the same method than the given fa.
#' @return The fa object, with in addition to the original list, the method used to compute the confidence interval, and 
#' a matrix containing estimated distribution parameters for each resample,
BootstrapCI <- function(fa, n.resamples=2.5E5, use_lmom = TRUE) {
  
  # extract fitted model parameters and flag as to whether the 
  # distribution is based on log transformed data
  base.params <- fa$distribution$parameters
  isTransformed <- fa$distribution$logTransformed
  # create output matrices to store parameter sets and quantile estimates
  if(fa$method == "mme") {nb_params = 3}
  else {nb_params = length(base.params$para)}
  param.sets <- matrix(NA, nrow = n.resamples, ncol = nb_params)
  if (fa$method == "lmom" || use_lmom)
  {
    ci_bootstrap = "lmom"

   
    # begin bootstrapping procedure
    for(i in 1:n.resamples) 
    {
      
      valid.moments <- FALSE
      j <- 0
      
      # allow up to 20 re-tries to re-sample 
      while(!valid.moments & j < 20) {  
        
        # sample 'n' random variates from base distribution
        if(!fa$method == "mme")
        {
          data <- rlmomco(n=length(fa$data), base.params)
        } else
        {
          data <- PearsonDS::rpearsonIII(length(fa$data),
                                         fa$distribution$parameters$estimate[1],
                                         fa$distribution$parameters$estimate[2],
                                         fa$distribution$parameters$estimate[3])
        }
        

        # compute sample l-moments
        sample.moms = lmom.ub(data)
        valid.moments <- are.lmom.valid(sample.moms)
      
        j <- j + 1
      }
      
      # error handling
      if(!valid.moments) {
        stop("Bootstrapping failed to sample valid moments")
      } else 
      {
        if(!fa$method == "mme")
        { type = base.params$type}
        else
        { type = "pe3"}
        # estimate distribution parameters
        dist.par <- lmom2par(sample.moms, type)
        
        # store the distribution parameters
        param.sets[i,] <- dist.par$para
      } 
    }  
    
  
  }
  else if (fa$method == "mle")
  {
    ci_bootstrap = "mle"
   
    # begin bootstrapping procedure
    for(i in 1:n.resamples) 
    {
      
      valid.moments <- FALSE
      j <- 0
      
      # allow up to 20 re-tries to re-sample 
      while(!valid.moments & j < 20) {  
        
        # sample 'n' random variates from base distribution
        data <- rlmomco(n=length(fa$data), base.params)
        init <- lmoms(data)
        init$ratios[3] <- 0 
        init <- fa$distribution$parameters
        distribution <-  fa$distribution$parameters$type
        dist.par <- mle2par(data, type=distribution, init.lmr=init) 
        valid.moments <- !is.null(dist.par)
        
        j <- j + 1
      }
      
      # error handling
      if(!valid.moments) {
        stop("Bootstrapping failed to sample valid moments")
      } else 
      {
        param.sets[i,] <- dist.par$para
      } 
    }  
      
  } else if (fa$method == "mme" && fa$distribution$name == "pe3")
  {
    ci_bootstrap = "mme"
    param.sets <- matrix(NA, nrow = n.resamples, ncol = 3)

    if (fa$distribution$logTransformed) 
    {
      series = log10(fa$data)
    } else
    {
      series = fa$data
    }
    
    
    ePIII <- function(x, order){
      # compute (centered) empirical centered moments of the data
      if(order == 1) return(mean(x))
      if(order == 2) return(stats::var(x))
      if(order == 3) return(e1071::skewness(x, type = 2))
    }
    # begin bootstrapping procedure
    for(i in 1:n.resamples) {
      
      valid.moments <- FALSE
      j <- 0
      
      # allow up to 20 re-tries to re-sample 
      while(!valid.moments & j < 20) {  
        
        # sample 'n' random variates from base distribution
        data <- PearsonDS::rpearsonIII(length(series),
                                       fa$distribution$parameters$estimate[1],
                                       fa$distribution$parameters$estimate[2],
                                       fa$distribution$parameters$estimate[3])
        
        m <- mean(data)
        v <- stats::var(data)
        s <- stats::sd(data)
        g <- e1071::skewness(data, type = 2)
        my_shape <- (2 / g) ^ 2
        my_scale <- sqrt(v) / sqrt(my_shape) * sign(g)
        my_location <- m - my_scale * my_shape
        
        start <- list(shape = my_shape, location = my_location, scale = my_scale)
        fit_mme <- fitdistrplus::fitdist(data, 'PIII', start = start, calcvcov = FALSE,
                                         method = "mme", order = 1:3, memp = ePIII, control = list(maxit = 1000))
        valid.moments <- !is.null(fit_mme)
        j <- j + 1
      }
      
      # error handling
      if(!valid.moments) {
        stop("Bootstrapping failed to sample valid moments")
      } else {
        param.sets[i,] <- fit_mme$estimate
  

        
        
      } 
      
    }
      
  }
  else
  {
    stop(
      sprintf('Method not \'%s\' not recognized!', distribution))
  }

  # now return list object containing output
  fa$ci.parameters <- param.sets
  fa$ci.bootstrap <- ci_bootstrap 

  return(fa)
}



#' FrequencyPlot Generates a nice-looking (hydrologist-centric) frequency plot
#' @param fa the result of FA_fit, with optionally the result from BootstrapCI
#' @export ggplot
#' @import ggplot2
#' @import scales
frequencyPlot <- function(fa) {
  
  library(ggplot2)
  library(scales)
  
  series <- fa$data
  
  fa_plot <- FA_getvalues(fa, nep=nonexceeds())
  # determine plotting positions
  bwpeaks <- data.frame(PROB = pp(series, sort = FALSE), FLOW = series)
  xbreaks <- c(0.002,0.01,0.1,0.25,0.5,0.8,0.9,0.95,0.975,0.99,0.995, 0.998)
  log.range <- log10(range(series, fa_plot[,-2:-1], na.rm = TRUE))
  lower <- 10^floor(log.range[1])
  upper <- 10^ceiling(log.range[2])
  cap <- lower
  ybreaks <- NULL
  while(cap < upper) {
    ybreaks <- c(ybreaks, seq(cap, cap*9, by = cap))
    cap <- cap * 10
  }
  
  if(!is.null(fa$ci.parameters))
  {
    # now plot
    ggplot(bwpeaks) + 
      geom_point(aes(x=PROB, y=FLOW)) + 
      theme_bw() + 
      scale_y_continuous(trans="log10", breaks=ybreaks, name="Discharge") +
      scale_x_continuous(trans=probability_trans(distribution="norm"),
                         breaks=xbreaks, labels=signif(prob2T(xbreaks), digits=4),
                         name="Return period [yrs]") +
      geom_line(data=fa_plot, aes(x=nep, y=estimate), color="red") +
      geom_line(data=fa_plot, aes(x=nep, y=lower), color="red", lty=2) +
      geom_line(data=fa_plot, aes(x=nep, y=upper), color="red", lty=2)
  }
  else
  {
    # now plot
    ggplot(bwpeaks) + 
      geom_point(aes(x=PROB, y=FLOW)) + 
      theme_bw() + 
      scale_y_continuous(trans="log10", breaks=ybreaks, name="Discharge") +
      scale_x_continuous(trans=probability_trans(distribution="norm"),
                         breaks=xbreaks, labels=signif(prob2T(xbreaks), digits=4),
                         name="Return period [yrs]") +
      geom_line(data=fa_plot, aes(x=nep, y=estimate), color="red")
  }

  
}

#' FA_comparisons computes the frequency analyses of some data with several statistical distribution,
#' and return a dataframe with the goodness-of-fit estimators
#' @param series A vector representing an extreme value series (e.g., annual maximum flood)
#' @param method The method used to fit the distribution. Either lmom" for the L-moments method (default value) 
#' or mle" for the maximum likelihood estimation. This function won't work with the moment matching estimation
#' @param list_distr A list of three-character name of the extreme value distribution 
#' @return A dataframe with the distribution name and the goodness-of-fit estimator : 
#' Chi-square test (providing a Chi-square score and a p-value)
#' the Akaike Information Criterion (AIC) and its variant the Bayesian Information Criterion (BIC),
#' the Kolmogorov-Smirnov test (providing the Kolmogorov-Smirnov score and a p-value) 
#' and the Anderson-Darling test (also providing a score and a p-value). 
#' Also provides the distribution ranking for all the estimators

FA_comparisons <- function(series, method = "lmom", list_distr = c("wei", "gev", "gum", "lp3", "lno", "nor", "exp", "gam", "gno"), nb_bins = NULL)
{

  chi_square <- c()
  chi_square_p <- c()
  AIC <- c()
  BIC <- c()
  KS_D <- c()
  KS_p <- c()
  AD_A2 <- c()
  AD_p <- c()
  
  # A first histogram is computed by not displayed, just to measure the maximum density 
  if (!is.null(nb_bins))
  {
    histogram <- hist(flow, breaks = nb_bins, prob = TRUE, plot=FALSE)
  } else
  {
    histogram <- hist(flow, prob = TRUE, plot=FALSE)
  }
  max_density <- max(histogram$density)
  curves_xvalues <- seq(from = min(histogram$breaks), to = max(histogram$breaks), by = (max(histogram$breaks) - min(histogram$breaks))/101)

  i = 0
  list_fa <- list()
  for (distr in list_distr)
  {
    # do the Frequency Analysis and store fitting metrics
    i <- i + 1
    fa <- FA_fit(series=series, distribution=distr, method=method)
    list_fa[[i]] <- fa
    chi_square <- c(chi_square, fa$chi_square$Chi_Square)
    chi_square_p <- c(chi_square_p, fa$chi_square$p_value)
    AIC <- c(AIC, fa$AIC)
    BIC <- c(BIC, fa$BIC)
    KS_D <- c(KS_D, fa$KS$statistic)
    KS_p <- c(KS_p, fa$KS$p.value)
    AD_A2 <- c(AD_A2, fa$AD$statistic)
    AD_p <- c(AD_p, fa$AD$p.value)
    
    # Check the maximum density for the fitted curve
    if (!fa$distribution$logTransformed)
    {
      fitted_curve <- dlmomco(curves_xvalues, fa$distribution$parameters)
    }
    else
    {
      # the pdf needs to be transformed 
      fitted_curve <- dlmomco(log10(curves_xvalues), fa$distribution$parameters)/(log(10)*curves_xvalues)
    }
    max_density <- max(max_density, max(fitted_curve))

  }

  # Let's display the histogram and the fitted curves
  # list of colors for the distributions
  colors = c(palette(), rainbow(8))
  # y axis is based on the maximum density on both histogram and curves
  ylim <- c(0, max_density*1.1)
  if (!is.null(nb_bins))
  {
    hist(flow, breaks = nb_bins, prob = TRUE, main = "Discharge Data with Fitted Distributions", xlab = "Discharge", ylim = ylim)
  } else
  {
    hist(flow, prob = TRUE, main = "Discharge Data with Fitted Distributions", xlab = "Discharge", ylim = ylim)
  }
  i <- 0
  for (distr in list_distr)
  {
    i <- i + 1
    fa <- list_fa[[i]]
    if (!fa$distribution$logTransformed)
    {
      curve(dlmomco(x, fa$distribution$parameters), add = TRUE, col = colors[i], lwd = 2)
    }
    else
    {
      # the pdf needs to be transformed 
      curve(dlmomco(log10(x), fa$distribution$parameters)/(log(10)*x), add = TRUE, col = colors[i], lwd = 2)
    }
  }
  # Legend
  legend("topright", legend = list_distr, col = colors, lwd = 2)
  
  results <- data.frame(
    distribution = list_distr,
    chi_square = chi_square,
    chi_square_rank = rank(chi_square),
    chi_square_p = chi_square_p,
    AIC = AIC,
    AIC_rank = rank(AIC),
    BIC = BIC,
    BIC_rank = rank(BIC),
    KS_D = KS_D,
    KS_p = KS_p,
    KS_rank = rank(KS_D),
    AD_A2 = AD_A2,
    AD_p = AD_p,
    AD_rank = rank(AD_A2)
    )
  
  return(results)
  
}

