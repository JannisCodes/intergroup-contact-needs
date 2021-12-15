#obj <- mdlWorker$lmerInterceptAttTypeTbl
#objz <- mdlWorker$lmerInterceptAttTypeTblZ

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
  confZ <- confint(objz, 
                   method = "boot", 
                   parallel = "multicore", 
                   ncpus = parallel::detectCores(logical = TRUE)) %>% 
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
        mutate(coef = gsub('_zwc|Z$', '', coef))
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
    tau00 = ifelse(length(tau00)>1, paste(tau00, collapse = "\n"), tau00),
    tau11 = ifelse(length(tau11)>1, paste(tau11, collapse = "\n"), tau11),
    rho01 = ifelse(length(rho01)>1, paste(rho01, collapse = "\n"), rho01),
    rho00 = ifelse(length(rho00)>1, paste(rho00, collapse = "\n"), rho00)
  )
  
  
  ## Model fit statistics.
  data = deparse(as.list(obj@call)$data)
  formula = Reduce(paste, deparse(as.list(obj@call)$formula))
  ll <- logLik(obj)[1]
  deviance <- deviance(obj)
  AIC <- AIC(obj)
  BIC <- BIC(obj)
  N <- as.numeric(obj@devcomp$dims["n"])
  G <- as.numeric(as.data.frame(smry$gvars)["# groups"])
  ICC <- performance::icc(obj)$ICC_adjusted
  rSqMarg <- as.numeric(performance::r2(obj)$R2_marginal) # Marginal R2
  rSqCond <- as.numeric(performance::r2(obj)$R2_conditional) # Conditional R2
  R2Marg_R2Cond <- paste(format(round(rSqMarg, 3), nsmall = 3), 
                         format(round(rSqCond, 3), nsmall = 3), 
                         sep = " / ")
  sumstat <- list(
    data = data,
    formula = formula,
    groupId = groupNam,
    logLik = ll, 
    deviance = deviance, 
    AIC = AIC,
    BIC = BIC, 
    N = N, 
    Groups = G,
    ICC = ICC,
    R2_marginal = rSqMarg,
    R2_conditional = rSqCond,
    R2Marg_R2Cond = R2Marg_R2Cond
  )
  
  # Combined Table output
  mdlCoefTbl <- data.frame(
    coef = coefOut$coef,
    B = coefOut$B,
    Beta = coefOut$Beta
  )
  mdlRandomTbl <- t(as.data.frame(restat)) %>%
    magrittr::set_colnames(c("B")) %>%
    as.data.frame %>%
    tibble::rownames_to_column(., "coef") %>%
    mutate(B = gsub("\\(\\)|\\(\\.\\)", "", B),
           Beta = "")
  mdlFitTbl <- t(as.data.frame(sumstat)) %>%
    magrittr::set_colnames(c("B")) %>%
    as.data.frame %>%
    tibble::rownames_to_column(., "coef") %>%
    mutate(Beta = "")
  mdlTbl <- rbind(mdlCoefTbl, mdlRandomTbl, mdlFitTbl)
  
  out <- list(
    coeficients = coefOut,
    random = restat,
    fit = sumstat,
    mdlTbl = mdlTbl
  )
  out
}
