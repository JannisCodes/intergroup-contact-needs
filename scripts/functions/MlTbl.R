#obj <- tblLmerS3AttSlopesAllportCore
#objz <- tblLmerS3AttSlopesAllportCoreZ

lmerTblPrep <- function (obj, objz = NULL, alpha = 0.05, ...) {
  # Summarize models
  smry <- summ(obj)
  smryZ <- summ(objz)
  
  # save grouing variable
  groupNam <- as.character(as.data.frame(smry$gvars)["Group"])
  
  # Get standardized coefficients 
  coefZ <- smryZ$coeftable %>% 
    as.data.frame %>%
    tibble::rownames_to_column(., "coef")
  confZ <- confint(objz, method = "boot") %>% 
    as.data.frame %>%
    tibble::rownames_to_column(., "coef")
  coefZOut <- merge(coefZ, confZ)
  colnames(coefZOut) <- c("coef", "est", "se", "tval", "df", "p", "lwr", "upr")
  coefZOut$Beta <- paste0(
    format(round(coefZOut$est, 2), nsmall = 2),
    " [", 
    format(round(coefZOut$lwr, 2), nsmall = 2),
    ", ",
    format(round(coefZOut$upr, 2), nsmall = 2),
    "]"
  )
  
  # Get unstandardized coefficients
  coef <- smry$coeftable %>% 
    as.data.frame %>%
    tibble::rownames_to_column(., "coef")
  conf <- confint(obj) %>% 
    as.data.frame %>%
    tibble::rownames_to_column(., "coef")
  coefOut <- merge(coef, conf)
  colnames(coefOut) <- c("coef", "est", "se", "tval", "df", "p", "lwr", "upr")
  coefOut$B <- paste0(
    format(round(coefOut$est, 2), nsmall = 2),
    ifelse(coefOut$p < .0001, "****", 
           ifelse(coefOut$p < .001, "***", 
                  ifelse(coefOut$p < .01, "**", 
                         ifelse(coefOut$p < .05, "*", "")))),
    " [", 
    format(round(coefOut$lwr, 2), nsmall = 2),
    ", ",
    format(round(coefOut$upr, 2), nsmall = 2),
    "]"
  )
  coefOut <- 
    merge(
      coefOut, 
      coefZOut %>% 
        select(coef, Beta) %>% 
        mutate(coef = gsub('Z$', '', coef))
    )
  
  
  ## Random Effects
  objVar <- insight::get_variance(obj)
  # Residual variance
  sigmaSq <- objVar$var.residual 
  # Random intercept variance
  tau00 <-  paste0(
    format(round(as.numeric(objVar$var.intercept), 2), nsmall = 2), 
    " (", names(objVar$var.intercept), ")"
  ) 
  # Random slope variance
  tau11 <- paste0(
    format(round(as.numeric(objVar$var.slope), 2), nsmall = 2), 
    " (", names(objVar$var.slope), ")"
  ) 
  # random-slope-intercept-correlation
  rho01 <- paste0(
    format(round(as.numeric(objVar$cor.slope_intercept), 2), nsmall = 2), 
    " (", names(as.data.frame(objVar$cor.slope_intercept)),".", row.names(objVar$cor.slope_intercept), ")"
  ) 
  # correlation between random slopes
  rho00 <- paste0(
    format(round(as.numeric(objVar$cor.slopes), 2), nsmall = 2), 
    " (", names(objVar$cor.slopes), ")"
  )  
  # collect random effect statistics
  restat <- list(
    sigmaSq = sigmaSq,
    tau00 = tau00,
    tau11 = tau11,
    rho01 = rho01,
    rho00 = rho00
  )
  
  
  ## Model fit statistics.
  ll <- logLik(obj)[1]
  deviance <- deviance(obj)
  AIC <- AIC(obj)
  BIC <- BIC(obj)
  N <- as.numeric(obj@devcomp$dims["n"])
  G <- as.numeric(as.data.frame(smry$gvars)["# groups"])
  ICC <- performance::icc(obj)$ICC_adjusted
  rSqMarg <- as.numeric(performance::r2(obj)$R2_marginal) # Marginal R2
  rSqCond <- as.numeric(performance::r2(obj)$R2_conditional) # Conditional R2
  sumstat <- list(
    groupId = groupNam,
    logLik = ll, 
    deviance = deviance, 
    AIC = AIC,
    BIC = BIC, 
    N = N, 
    Groups = G,
    ICC = ICC,
    R2_marginal = rSqMarg,
    R2_conditional = rSqCond
  )
  
  out <- list(
    coeficients = coefOut,
    random = restat,
    fit = sumstat
  )
  out
}
