## Author: Jack Pohlmann
## Date: July 19th, 2023


#### Package Dependencies####
# library(tidyverse)
# library(magrittr)
# library(lubridate)
# library(readxl)
# library(kableExtra)
# library(nnet) #multinomial log-reg model


#### Table Functions ####
p <- function(beta, stdev){
  ## returns the p-value given a beta and stdev
  zscore <- beta/stdev
  pval <- (1 - pnorm(abs(zscore), 0, 1))*2 
  return(pval)
}

pval_format <- function(pval){
  ## Formats p-values to appropriate decimal places
  if (is.na(pval)){
    return('NA')
  } else if (pval<0.01){
    return('<0.01')
  }
  else{
    return(as.character(round(pval, 2)))
  }
}

beta_format <- function(beta, digits=2){
  ## Formats beta values to appropriate decimal places
  beta <- round(beta, digits = digits)
  return(beta)
}

adj_pct <- function(x,y){
  ## Function to change the % of a variable
  pct <- round(x/y*100,1)
  out <- paste0(x, ' (', pct, '%)')
  return(out)
}

combined_tbl <- function(y, vars, data, stat='mean'){
  #-----------------------------------------------------------------------------
  # Function designed to create semi-customized Table 1s, with counts for both  
  # total cohort and stratified groups.
  #
  # y: str, name of column in the supplied data frame by which you would like to  
  #    stratify. Note that we should have type(data$y) == factor
  #
  # vars: vector of column names for the covariates of interest. Columns for the  
  #       co-variates should be either numeric, integer, or factor.
  #
  # data: data frame containing columns y and vars.
  #
  # NOTE: This function will spit out a raw, ugly table. The user will need to manually edit
  # the table after the function is run to put it in a more presentable format. I recommend
  # using the "kable" and "kableExtra" libraries to knit beautiful HTML and PDF tables. See
  # their respective documentation for more information.
  # 
  #-----------------------------------------------------------------------------
  
  d <- y
  y <- data[[y]]
  o_tbl <- table(y) ## get counts of patients in each stratum
  N <- nrow(data)
  
  ## set up empty table, add total row
  df <- data.frame(matrix(ncol = 3+length(levels(y)), nrow = 0))
  out.row <- c('Total', as.character(N), paste0(o_tbl, ' (', round(o_tbl/N,3)*100, '%)'), '-')
  df <- rbind(df, out.row)
  
  if (stat=='mean'){

  ## iterate through the list of variables provided
  for (i in 1:length(vars)) {
    
    if( !vars[i] %in% colnames(data) ){
      print(paste0('The variable ',vars[i],' was not found in the dataframe.'))
    } 
    else {
      x <- data[[vars[i]]]
      
      ## numeric and integer variables
      if (class(x)=='numeric' | class(x)=='integer') {
        u <- round(mean(x, na.rm = T), 1)
        stdev <- round(sd(x, na.rm = T), 1)
        out.row <- c(vars[i], paste0(u,' (', stdev,')'))
        for (j in 1:length(levels(y))) {
          lvl <- levels(y)[j]
          u <- round(mean(x[which(y==lvl)], na.rm = T), 1)
          stdev <- round(sd(x[which(y==lvl)], na.rm = T), 1)
          out <- paste0(u, ' (', stdev, ')')
          out.row <- c(out.row, out)
        }
        
        ## determine the appropriate statistical test to use
        if (length(levels(y)) == 2){
          ## t-test
          pval <- pval_format(t.test(x[ which(y==levels(y)[1]) ], x[ which(y==levels(y)[2]) ])$p.value)
          out.row <- c(out.row, pval)
        } else if (length(levels(y)) > 2) { 
          ## anova
          pval <- pval_format(summary(aov(x~y))[[1]]$'Pr(>F)'[1])
          out.row <- c(out.row, pval)
        }
        df <- rbind(df, out.row)
        

      } 
      else if (class(x)=='factor') {
        n_tbl <- table(x, y)
        pval <- pval_format(chisq.test(n_tbl)$p.value)
        
        if ( identical(levels(x), c('0','1')) ){
            n <- length(which(x=='1'))
            pct_tbl <- round((n_tbl[2,]/o_tbl)*100, 1)
            out.row <- c(vars[i], paste0(n, ' (',round(n/N, 3)*100, '%)'),
                         paste0(n_tbl[2,], ' (', pct_tbl, '%)'),
                         pval)
            df <- rbind(df, out.row)
            }
          
        else {
          
          totVal <- length(which(!is.na(x)))
          ndata <- subset(data, !is.na(x)) 
          t_tbl <- table(ndata[[d]])
          
          if(totVal < N) {
            out.row <- c(vars[i], paste0(totVal, ' (',round(totVal/N, 3)*100, '%)'),
                         paste0(t_tbl, ' (', round(t_tbl/o_tbl,3)*100, '%)'), 
                         pval)
          }else{
            out.row <- c(vars[i],
                         rep('', 1 + length(levels(y))), pval)
          }
          df <- rbind(df, out.row)
          
          for (j in 1:length(levels(x))) {
            lvl <- levels(x)[j]
            n <- length(which(x==lvl))
            pct_tbl <- round((n_tbl[j,]/t_tbl)*100, 1)
            out.row <- c(paste0( levels(x)[j]), paste0(n, ' (',round(n/totVal, 3)*100, '%)'),
                         paste0(n_tbl[j,], ' (', pct_tbl, '%)'),
                         '')
            df <- rbind(df, out.row)
          }}}
      else {
        print(paste0('Variable ',vars[i], 'is neither numeric or a factor.'))
      }
    }
  }
  
  ## rename the table columns
  names(df) <- c('Variable', 'Count(%)/Mean(SD)',levels(y), 'p-value')
  
  } 
  else if (stat=='median'){
    
    ## iterate through the list of variables provided
    for (i in 1:length(vars)) {
      x <- data[[vars[i]]]
      
      ## dummy line for headers
      if(vars[i] == "")
      {
        df <- rbind(df, c('temp', rep('', 5)))
      }
      else {
        ## numeric and integer variables
        if (class(x)=='numeric' | class(x)=='integer') {
          qs <- quantile(x, c(.25, .50, .75), na.rm = T)
          m <- round(qs[2],2)
          q1 <- round(qs[1], 2)
          q3 <- round(qs[3], 2)
          out.row <- c(vars[i], paste0(m,' (',q1,'-',q3,')'))
          
          for (j in 1:length(levels(y))) {
            lvl <- levels(y)[j]
            qs <- quantile(x[which(y==lvl)], c(.25, .50, .75), na.rm = T)
            m <- round(qs[2],2)
            q1 <- round(qs[1], 2)
            q3 <- round(qs[3], 2)
            out <- c(paste0(m,' (',q1,'-',q3,')'))
            
            out.row <- c(out.row, out)
          }
          
          ## determine the appropriate statistical test to use
          # To do...
          
          df <- rbind(df, out.row)
          
          ## factor variables
        } else if (class(x)=='factor') {
          n_tbl <- table(x, y)
          pval <- pval_format(chisq.test(n_tbl)$p.value)
          
          if ( identical(levels(x), c('0','1')) ){
                n <- length(which(x=='1'))
                pct_tbl <- round((n_tbl[2,]/o_tbl)*100, 1)
                out.row <- c(vars[i], paste0(n, ' (',round(n/N, 3)*100, '%)'),
                             paste0(n_tbl[2,], ' (', pct_tbl, '%)'),
                             pval)
                df <- rbind(df, out.row)
            }
          else {
            
            totVal <- length(which(!is.na(x)))
            ndata <- subset(data, !is.na(x)) 
            t_tbl <- table(ndata[[d]])
            
            if(totVal < N) {
              out.row <- c(vars[i], paste0(totVal, ' (',round(totVal/N, 3)*100, '%)'),
                           paste0(t_tbl, ' (', round(t_tbl/o_tbl,3)*100, '%)'), 
                           pval)
            }else{
              out.row <- c(vars[i],
                           rep('', 1 + length(levels(y))), pval)
            }
            df <- rbind(df, out.row)
            
            for (j in 1:length(levels(x))) {
              lvl <- levels(x)[j]
              n <- length(which(x==lvl))
              pct_tbl <- round((n_tbl[j,]/t_tbl)*100, 1)
              out.row <- c(paste0( levels(x)[j]), paste0(n, ' (',round(n/totVal, 3)*100, '%)'),
                           paste0(n_tbl[j,], ' (', pct_tbl, '%)'),
                           '')
              df <- rbind(df, out.row)
            }}}}}
    
    ## rename the table columns
    names(df) <- c('Variable', 'Count(%)/Median(Q1-Q3)',levels(y))
  } 
  else {
    print("Argument provided for 'stat' not recognised. Please choose either 'mean' or 'median'.")
  }
  return(df)
}

univar_numeric <- function(df, y, xcols){
  
  ## Create empty table to populate
  tab.out <- data.frame(matrix(ncol = 3*( length(levels(df[[y]]))-1 )+1, nrow = 0))
  
  ## iterate through the list of variables
  for (i in 1:length(xcols)){
    x <- xcols[i]
    
    #dummy line for headers
    if(x == "")
    {
      tab.out <- rbind(tab.out, c(rep('', 3*(length(levels(df[[y]]))-1)+1)))
    }
    else
    {
      ## run a regression with the outcome and the predictor
      mlog.out <- multinom(df[[y]] ~ df[[x]], trace=F)
      sum.out <- summary(mlog.out)
      beta <- sum.out$coefficients
      stdev <- sum.out$standard.errors
      # odds <- exp(beta)
      # ci.025 <- exp(beta-stdev)
      # ci.975 <- exp(beta+stdev)
    
      ## iterate through the levels of the outcome
      row.out <- c(x)
      if (length(sum.out$lev) == 2) #there are two levels
      {
        bet <- beta[2]
        std <- stdev[2]
        pval <- p(bet, std)
        or <- exp(bet)
        ci.025 <- exp(bet-1.96*std)
        ci.975 <- exp(bet+1.96*std)
        
        #format how values are displayed
        or <- round( or, 2 )
        if (pval < 0.05) {
          or <- paste0(or, '*')}
        ci.025 <- round(ci.025, 2)
        ci.975 <- round(ci.975, 2)

        row.out <- c(row.out,
                     bet,
                     std,
                     pval)
      }
      
      else #more than two levels
      {
      for (j in 1:(length(sum.out$lev )-1))
        {
        bet <- beta[j,2]
        std <- stdev[j,2]
        pval <- p(bet, std)
        or <- exp(bet)
        ci.025 <- exp(bet-1.96*std)
        ci.975 <- exp(bet+1.96*std)
        
        #format how values are displayed
        or <- round( or, 2 )
        if (pval < 0.05) {
          or <- paste0(or, '*')}
        ci.025 <- round(ci.025, 2)
        ci.975 <- round(ci.975, 2)
        
        row.out <- c(row.out,
                     or,
                     ci.025,
                     ci.975)
        }
      }
      tab.out <- rbind(tab.out, row.out)
    }
  }
  
  ## rename the table columns
  names(tab.out) <- c( 'Variable', rep(c('OR','CI 2.50%', 'CI 97.5%'), length(levels(df[[y]]))-1 ))
  
  return(tab.out)
}

univar_factor <- function(df, y, xcols){
  
  ## create an empty table
  tab.out <- data.frame(matrix(ncol = 3*( length(levels(df[[y]]))-1 )+1, nrow = 0))
  
  ## iterate through the list of predictor variables
  for (i in 1:length(xcols)){
    x <- xcols[i]
    
    if ( x == '' ){
      tab.out[nrow(tab.out) + 1, ] <- c(x, rep('', 6))
    } else
    {
    ## if the supplied variable isnt a factor type, make it one
    if (class(df[[x]]) != 'factor'){
      df[[x]] <- factor(df[[x]])}
    
    lvls <- levels(df[[x]])
    
    if( ((length(lvls) != 2 | x == 'ethnicity') & 
        length(which(is.na(df[[x]]))) == 0) | x == 'tici' | 
        x == 'lvo_location'){
      tab.out[nrow(tab.out) + 1, ] <- c(x, rep('', 6))
    }
    
    ## run a regression with the outcome y with the predictor x
    mlog.out <- multinom(df[[y]] ~ df[[x]], trace=F)
    
    
    sum.out <- summary(mlog.out)
    bet <- sum.out$coefficients
    std <- sum.out$standard.errors
    or <- exp(coef(mlog.out))
    ci <- exp(confint(mlog.out))
    
    ## iterate through the levels of the outcome
    for (j in 2:length(lvls)){
      r <- c(lvls[j])
      ## iterate through the levels of the predictor
      for ( k in 1:(length( sum.out$lev )-1) ){
        pval <- p(bet[k,j], std[k,j])
        or_format <- round(or[k,j],2)
        if (pval < 0.05){
          or_format <- paste0(or_format, '*')}
        
        r <- c(r, or_format, round(ci[j,1,k],2), round(ci[j,2,k],2))}
      tab.out[nrow(tab.out) + 1,] <- r}}
  }
  
  ## rename the table columns
  names(tab.out) <- c( 'Variable', rep(c('OR','CI 2.50%', 'CI 97.5%'), length( sum.out$lev )-1) )
  
  return(tab.out)
}

univar_analysis <- function(df, y, xcols){
  ## create an empty table
  y_name <- y
  y <- df[[y_name]]
  width <- if( is.null(levels(y)) ){3 +1} else {3*( length(levels(y))-1 )+1}
  tab.out <- data.frame(matrix(ncol = width, nrow = 0))
  
  ## iterate through the list of predictor variables
  for (i in 1:length(xcols)){
    
    #**********************************************
    break('IN DEVELOPMENT. FUNCTION INCOMPLETE.')
    #**********************************************
    
    x_name <- xcols[i]
    x <- df[[x_name]]
    
      ## if the supplied variable isnt a factor type, make it one
      if (class(x) == 'character' | class(x) == 'factor'){
        x <- as.factor(x)
        x_lvls <- levels(x)
      }
        
        tab.out[nrow(tab.out) + 1, ] <- c( x_name, rep('', width-1) )
        
        ## run a regression with the outcome y with the predictor x
        if( class(y) %in% c('numeric','integer') & !identical(unique(y), c(0,1)) ){ 
          #linear regression
          mod.out <- lm(y ~ x)
        } 
        else if ( length(levels(y)) == 2 | identical(unique(y), c(0,1)) ){
          ## log reg
          mod.out <- glm(y ~ x, family = binomial())
        }
        else if ( levels(y) > 2 ){
          ## multinom log reg
          mod.out <- multinom(y ~ x, trace=F)
        }
        
        sum.out <- summary(mod.out)
        bet <- sum.out$coefficients
        std <- sum.out$standard.errors
        or <- exp(coef(mod.out))
        ci <- exp(confint(mod.out))
        
        ## iterate through the levels of the outcome
        for (j in 2:length(lvls)){
          r <- c(lvls[j])
          ## iterate through the levels of the predictor
          for ( k in 1:(length( sum.out$lev )-1) ){
            pval <- p(bet[k,j], std[k,j])
            or_format <- round(or[k,j],2)
            if (pval < 0.05){
              or_format <- paste0(or_format, '*')}
            
            r <- c(r, or_format, round(ci[j,1,k],2), round(ci[j,2,k],2))}
          tab.out[nrow(tab.out) + 1,] <- r}
  
      
      
  }
  
  ## rename the table columns
  names(tab.out) <- c( 'Variable', rep(c('OR','CI 2.50%', 'CI 97.5%'), length( sum.out$lev )-1) )
  
  return(tab.out)
}

mvlr_tbl <- function(sumTbl){
  dat <- sumTbl[-1, c('term', 'estimate', 'std.error', 'p.value')]
  
  dat$Variable <- dat$term
  dat$`Odds Ratio` <- round(exp(dat$estimate),2)
  
  z <- 1.96
  ci_lo <- round(exp(dat$estimate - z*dat$std.error),2)
  ci_hi <- round(exp(dat$estimate + z*dat$std.error),2)
  ci <- paste0(ci_lo, ' - ', ci_hi)
  dat$`95% Conf Int` <- ci
  
  dat$`p-value` <- unlist(lapply(dat$p.value, pval_format))
  
  out <- dat[,c('Variable', 'Odds Ratio', '95% Conf Int', 'p-value')]
  row.names(out) <- NULL
  
  return(out)
}
