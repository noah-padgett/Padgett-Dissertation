# easily grab data for all pages

get_date <- function(){
  format(Sys.time(), "%Y_%m_%d")
}


# The follow check the normality and variance assumptions of ANOVA
# returns the summary stats by cell
anova_assumptions_check <- function(
  dat, outcome, factors, model = NULL, stats = NULL)
{
  if(is.null(stats) == T){
    g <- paste0('group',1:length(factors))
    stats <- c(g, 'n', 'mean', 'sd', 'skew',
               'kurtosis', 'min', 'Q0.25','median',
               'Q0.75','max')
  }
  ## Check to see if a model was supplied
  if(is.null(model) == T){
    model <- as.formula(paste(outcome, '~', paste(factors,collapse = "*")))
  }

  cat('\n ============================= \n')
  cat('\n Tests and Plots of Normality:\n')
  # Assess normality
  aov.out = aov(model, data = dat)
  plot(aov.out,ask=F)
  # shapiro-wilks test
  if(length(aov.out$residuals) > 5000){
    res <- sample(aov.out$residuals, 5000)
  } else res <- aov.out$residuals
  cat('\n Shapiro-Wilks Test of Normality of Residuals:\n')
  print(shapiro.test(res))

  # K-S Test
  cat('\n K-S Test for Normality of Residuals:\n')
  print(ks.test(aov.out$residuals, 'pnorm',
                alternative = 'two.sided'))
  cat('\n')

  # Histograms
  ## loop around ggplot to make histograms
  for(i in 1:length(factors)){
    print(
      ggplot(dat, aes_string(x=outcome,
                             fill = factors[i], color = factors[i])) +
        geom_histogram(alpha=0.5, position="identity") +
        labs(title = paste(outcome, 'distribution by', factors[i])) +
        theme_classic()
    ) ## Print plot
  }

  cat('\n ============================= \n')
  cat('\n Tests of Homogeneity of Variance\n')
  # Varainces
  ## loop around levene's test

  for(i in 1:length(factors)){
    cat('\n \n Levenes Test: ', factors[i], '\n \n \n')
    out <- leveneTest(as.formula(paste(outcome, '~',factors[i])),
                     data = dat,
                     center="mean")
    print(out)
  }
  # return(summary_stats)
}


## General Form of Omega
o2 <- function(ss_num,ss_dem, df, n, mse)
{
  (ss_num - df*mse)/(ss_dem + (n-df)*mse)
}
# Estimating Omega**2
omega2 <- function(fit.sum){
  k <- length(fit.sum[[1]]$`Mean Sq`)
  ms_error <- fit.sum[[1]]$`Mean Sq`[k]
  omega <- matrix(NA,ncol=1,nrow=(k-1))
  dfs <- fit.sum[[1]]$Df
  ss <- fit.sum[[1]]$`Sum Sq`
  sst <- sum(fit.sum[[1]]$`Sum Sq`)
  for(i in 1:(k-1))
  {
    omega[i,] <- round(o2(ss[i],sst, dfs[i], dfs[i]+1, ms_error),4)
  }
  rownames(omega) <- rownames(fit.sum[[1]])[1:(k-1)]
  colnames(omega) <- "omega^2"
  return(omega)
}
# Estimate partial-omega**2
p_omega2 <- function(fit.sum){
  k <- length(fit.sum[[1]]$`Mean Sq`)
  ms_error <- fit.sum[[1]]$`Mean Sq`[k]
  N <- sum(fit.sum[[1]]$Df) + 1
  p.omega <- matrix(NA,ncol=1,nrow=(k-1))
  dfs <- fit.sum[[1]]$Df
  sss <- fit.sum[[1]]$`Sum Sq`
  for(i in 1:(k-1))
  {
    p.omega[i,] <- round(o2(sss[i],sss[i], dfs[i], N, ms_error),4)
  }
  rownames(p.omega) <- rownames(fit.sum[[1]])[1:(k-1)]
  colnames(p.omega) <- "partial-omega^2"
  return(p.omega)
}


# custom boxplot
cust.boxplot <- function(x) {
  r <- quantile(x, probs = c(0.10, 0.25, 0.5, 0.75, 0.90))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}
