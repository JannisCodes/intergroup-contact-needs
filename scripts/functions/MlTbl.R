library(plyr)
#obj <- tblLmerS1AttSlopeCoreQlt
#objz <- tblLmerS1AttSlopeCoreQltZ

# obj <- mdlTblElements[[i]][[1]]
# objz <- mdlTblElements[[i]][[2]]

lmerTblPrep <- function (obj, objz = NULL, name = NULL, alpha = 0.05, ...) {
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
                   method = "Wald") %>% 
                   #method = "boot", 
                   #parallel = "multicore", 
                   #ncpus = parallel::detectCores(logical = TRUE)) %>% 
    as.data.frame %>%
    tibble::rownames_to_column(., "coef")
  coefZCI <- merge(coefZ, confZ)
  colnames(coefZCI) <- c("coef", "est", "se", "tval", "df", "p", "lwr", "upr")
  
  coefZR2 <- r2glmm::r2beta(
    model = objz,
    partial = TRUE,
    # 'sgv' standardized generalized variance approach, 
    # 'kr' (Kenward Roger approach; only with lme),
    # 'nsj' (Nakagawa and Schielzeth approach)
    method = 'sgv' 
  ) %>%
    as.data.frame %>%
    select(
      coef = Effect,
      Rsq,
      RsqLower = lower.CL,
      RsqUpper = upper.CL
    )
  coefZOut <- merge(coefZCI, coefZR2)
  coefZOut$Beta <- paste0(
    format(round(coefZOut$est, 2), nsmall = 2),
    " [", 
    format(round(coefZOut$lwr, 2), nsmall = 2),
    ", ",
    format(round(coefZOut$upr, 2), nsmall = 2),
    "]"
  )
  coefZOut$BetaStar <- ifelse(coefZOut$p < .0001, "****", 
                          ifelse(coefZOut$p < .001, "***", 
                                 ifelse(coefZOut$p < .01, "**", 
                                        ifelse(coefZOut$p < .05, "*", ""))))
  
  # Get unstandardized coefficients
  coef <- smry$coeftable %>% 
    as.data.frame %>%
    tibble::rownames_to_column(., "coef")
  conf <- confint(obj,
                  method = "Wald") %>% 
                  #method = "profile") %>% 
    as.data.frame %>%
    tibble::rownames_to_column(., "coef")
  coefOut <- plyr::join(coef, conf, by = "coef")
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
  coefOut$Bstar <- ifelse(coefOut$p < .0001, "****", 
                          ifelse(coefOut$p < .001, "***", 
                                 ifelse(coefOut$p < .01, "**", 
                                        ifelse(coefOut$p < .05, "*", ""))))
  coefOut <- 
    plyr::join(
      coefOut %>% 
        mutate(coef = gsub('_cwc|C$', '', coef)) %>%
        mutate(coef = gsub('C:', ':', coef)), 
      coefZOut %>% 
        select(coef, estBeta = est, Beta, BetaStar) %>% 
        mutate(coef = gsub('_zwc|_gmz|Z$', '', coef)) %>%
        mutate(coef = gsub('Z:', ':', coef)) %>%
        mutate(Beta = ifelse(coef == "(Intercept)", "", Beta)),
      by = "coef"
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
  G <- as.character(as.data.frame(smry$gvars)["# groups"])
  ICC <- ifelse(
    berryFunctions::is.error(performance::icc(obj)$ICC_adjusted), 
    performance::icc(obj),
    performance::icc(obj)$ICC_adjusted
  ) 
  rSqMarg <- as.numeric(performance::r2(obj)$R2_marginal) # Marginal R2
  rSqCond <- as.numeric(performance::r2(obj)$R2_conditional) # Conditional R2
  R2Marg_R2Cond <- paste(format(round(rSqMarg, 3), nsmall = 3), 
                         format(round(rSqCond, 3), nsmall = 3), 
                         sep = " / ")
  sumstat <- list(
    name = ifelse(is.null(name), "", name),
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
    beta = coefZOut,
    random = restat,
    fit = sumstat,
    mdlTbl = mdlTbl
  )
  out
}

#obj <- lmAllAttFreqQualX2
lmTblPrep <- function (obj, alpha = 0.05, ...) {
  # Summarize models
  smry <- summ(obj, confint = TRUE)
  smrySE <- summ(obj, confint = FALSE)
  
  # Get unstandardized coefficients
  coefCI <- smry$coeftable %>% 
    as.data.frame %>%
    tibble::rownames_to_column(., "coef")
  colnames(coefCI) <- c("coef", "est", "lwr", "upr", "tval", "p")
  
  coefSE <- smrySE$coeftable %>% 
    as.data.frame %>%
    tibble::rownames_to_column(., "coef") %>%
    select(coef, se = S.E.)
  
  coef <- plyr::join(
    coefCI, 
    coefSE,
    by = "coef"
  )
  
  coef$B <- paste0(
    format(round(coef$est, 2), nsmall = 2),
    ifelse(coef$p < .0001, "****", 
           ifelse(coef$p < .001, "***", 
                  ifelse(coef$p < .01, "**", 
                         ifelse(coef$p < .05, "*", "")))),
    " [", 
    format(round(coef$lwr, 2), nsmall = 2),
    ", ",
    format(round(coef$upr, 2), nsmall = 2),
    "]"
  )
  
  # calc effect size
  eta <- effectsize::eta_squared(obj, partial = TRUE) %>%
    as.data.frame %>%
    mutate(coef = Parameter) %>%
    mutate(Eta2_partial = format(round(Eta2_partial, 2), nsmall = 2))
  
  Rsq <- r2glmm::r2beta(
    model = obj,
    partial = TRUE,
    method = 'lm'
  ) %>% 
    as.data.frame %>%
    select(
      coef = Effect,
      Rsq,
      RsqLower = lower.CL,
      RsqUpper = upper.CL
    )
  
  # Merge B and effect size
  coefEta <- 
    plyr::join(
      coef %>% 
        mutate(coef = gsub('_c', '', coef)), 
       eta %>% 
        select(coef, Eta2_partial) %>% 
        mutate(coef = gsub('_c', '', coef)),
      by = "coef"
    )
  coefOut <- 
    plyr::join(
      coefEta, 
      Rsq %>% 
        mutate(coef = gsub('_c', '', coef)),
      by = "coef"
    )
  
  
  ## Model fit statistics.
  data = deparse(as.list(smry$model$call)$data)
  formula = Reduce(paste, deparse(as.list(smry$model$call)$formula))
  N <- nrow(smry$model$model)
  fTest <-
    paste(
      "F(",
      summary(obj)$fstatistic[2],
      ", ",
      summary(obj)$fstatistic[3],
      ") = ",
      format(round(summary(obj)$fstatistic[1], 2), nsmall = 2),
      ", ",
      ifelse(
        pf(
          summary(obj)$fstatistic[1],
          summary(obj)$fstatistic[2],
          summary(obj)$fstatistic[3],
          lower.tail = FALSE
        ) < .001,
        "p < .001",
        paste("p = ", format(round(
          pf(
            summary(obj)$fstatistic[1],
            summary(obj)$fstatistic[2],
            summary(obj)$fstatistic[3],
            lower.tail = FALSE
          ),
          3
        ), nsmall = 3))
      ),
      sep = ""
    )
  
  R2 <- as.numeric(performance::r2(obj)$R2) # Marginal R2
  R2_adjusted <- as.numeric(performance::r2(obj)$R2_adjusted) # Conditional R2
  R2_R2_adjusted <- paste(format(round(R2, 3), nsmall = 3), 
                         format(round(R2_adjusted, 3), nsmall = 3), 
                         sep = " / ")
  sumstat <- list(
    data = data,
    formula = formula, 
    N = N, 
    fTest = fTest,
    R2 = R2,
    R2_adjusted = R2_adjusted,
    R2_R2_adjusted = R2_R2_adjusted
  )
  
  # Combined Table output
  mdlCoefTbl <- data.frame(
    coef = coefOut$coef,
    B = coefOut$B,
    Eta2_partial = coefOut$Eta2_partial
  )
  mdlFitTbl <- t(as.data.frame(sumstat)) %>%
    magrittr::set_colnames(c("B")) %>%
    as.data.frame %>%
    tibble::rownames_to_column(., "coef") %>%
    mutate(Eta2_partial = "")
  mdlTbl <- rbind(mdlCoefTbl, mdlFitTbl)
  
  out <- list(
    coeficients = coefOut,
    fit = sumstat,
    mdlTbl = mdlTbl
  )
  out
}

