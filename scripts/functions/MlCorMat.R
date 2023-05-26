MlCorMat <-
  function(data,
           id,
           selection,
           labels = selection,
           method = c("pearson", "spearman"),
           removeTriangle = c("upper", "lower"),
           result = c("none", "html", "latex")) {
    # for testing
    # data <- dtWorkerSupp$workerInteractionType
    # id <- "PID"
    # selection <- varSelection
    # labels <- varlabels
    
    # calculate trait scores (participant means)
    data <- 
      data %>%
      select(id, selection) %>%
      group_by_at(vars(matches(id))) %>%
      mutate_at(selection, list(trait = ~mean(., na.rm=TRUE))) %>%
      ungroup 
    
    data <-
      map_dfc(selection,
              ~ data %>%
                transmute(
                  !!str_c(.x, '_state') :=
                    !!rlang::sym(.x)-!!rlang::sym(str_c(.x, "_trait"))
                )) %>%
      bind_cols(data, .)
    
    # Fid out why own is different
    # raw <-
    #   data %>%
    #   select(selection) %>%
    #   as.matrix
    # 
    # corRaw <- rcorr(raw, type=method[1])
    # rRaw <- corRaw$r # Matrix of correlation coeficients
    # pRaw <- corRaw$P # Matrix of p-value 
    # 
    # ## Define notions for significance levels; spacing is important.
    # starsRaw <- ifelse(pRaw < .001, "***", ifelse(pRaw < .01, "**", ifelse(pRaw < .05, "* ", "   ")))
    # 
    # ## trunctuate the correlation matrix to two decimal
    # rRaw <- format(round(cbind(rep(-1.11, ncol(raw)), rRaw), 2))[,-1]
    # 
    # ## build a new matrix that includes the correlations with their apropriate stars
    # rRawStar <- matrix(paste(rRaw, starsRaw, sep=""), ncol=ncol(raw))
    # 
    # 
    # trait <-
    #   data %>%
    #   select(paste0(selection, "_trait")) %>%
    #   as.matrix
    # 
    # corTrait <- rcorr(trait, type=method[1])
    # rTrait <- corTrait$r # Matrix of correlation coeficients
    # pTrait <- corTrait$P # Matrix of p-value 
    # 
    # ## Define notions for significance levels; spacing is important.
    # starsTrait <- ifelse(pTrait < .001, "***", ifelse(pTrait < .01, "**", ifelse(pTrait < .05, "* ", "   ")))
    # 
    # ## trunctuate the correlation matrix to two decimal
    # rTrait <- format(round(cbind(rep(-1.11, ncol(trait)), rTrait), 2))[,-1]
    # 
    # ## build a new matrix that includes the correlations with their apropriate stars
    # rTraitStar <- matrix(paste(rTrait, starsTrait, sep=""), ncol=ncol(trait))
    # 
    # 
    # state <-
    #   data %>%
    #   select(paste0(selection, "_state")) %>%
    #   as.matrix
    # 
    # corState <- rcorr(state, type=method[1])
    # rState <- corState$r # Matrix of correlation coeficients
    # pState <- corState$P # Matrix of p-value 
    # 
    # ## Define notions for significance levels; spacing is important.
    # starsState <- ifelse(pState < .001, "***", ifelse(pState < .01, "**", ifelse(pState < .05, "* ", "   ")))
    # 
    # ## trunctuate the correlation matrix to two decimal
    # rState <- format(round(cbind(rep(-1.11, ncol(state)), rState), 2))[,-1]
    # 
    # ## build a new matrix that includes the correlations with their apropriate stars
    # rStateStar <- matrix(paste(rState, starsState, sep=""), ncol=ncol(state))
    # 
    # ## combine state and trait
    # rComb <- matrix("", nrow = length(selection), ncol = length(selection))
    # rComb[upper.tri(rComb, diag = FALSE)] <- rTraitStar[upper.tri(rTraitStar, diag = FALSE)]
    # rComb[lower.tri(rComb, diag = FALSE)] <- rStateStar[upper.tri(rStateStar, diag = FALSE)]
    # rownames(rComb) <- selection
    # colnames(rComb) <- paste(selection, "", sep="")
    
    # Multilevel correlations based on ml model
    corML <- misty::multilevel.cor(data %>% select(selection), data %>% select(id), output = FALSE)
    
    # Within ppt
    rMlWit <- corML$result$with.cor
    pMlWit <- corML$result$with.p
    
    ## Define notions for significance levels; spacing is important.
    starMlWit <- ifelse(pMlWit < .001, "***", ifelse(pMlWit < .01, "**", ifelse(pMlWit < .05, "* ", "   ")))
    
    ## trunctuate the correlation matrix to two decimal
    rMlWit <- format(round(cbind(rep(-1.11, length(selection)), rMlWit), 2))[,-1]
    
    ## build a new matrix that includes the correlations with their apropriate stars
    rMlWitStar <- matrix(paste(rMlWit, starMlWit, sep=""), ncol=length(selection))
    
    # between ppt
    rMlBtw <- corML$result$betw.cor
    pMlBtw <- corML$result$betw.p
    
    ## Define notions for significance levels; spacing is important.
    starMlBtw <- ifelse(pMlBtw < .001, "***", ifelse(pMlBtw < .01, "**", ifelse(pMlBtw < .05, "* ", "   ")))
    
    ## trunctuate the correlation matrix to two decimal
    rMlBtw <- format(round(cbind(rep(-1.11, length(selection)), rMlBtw), 2))[,-1]
    
    ## build a new matrix that includes the correlations with their apropriate stars
    rMlWitBtw <- matrix(paste(rMlBtw, starMlBtw, sep=""), ncol=length(selection))
    
    rMlComb <- matrix("", nrow = length(selection), ncol = length(selection))
    rMlComb[upper.tri(rMlComb, diag = FALSE)] <- rMlWitBtw[upper.tri(rMlWitBtw, diag = FALSE)]
    rMlComb[lower.tri(rMlComb, diag = FALSE)] <- rMlWitStar[upper.tri(rMlWitStar, diag = FALSE)]
    rownames(rMlComb) <- labels
    colnames(rMlComb) <- paste(labels, "", sep="")
    
    ## Descriptives
    descriptives <- data.frame(
      variable = labels,
      `Grand Mean` = data %>%
        select(paste0(selection, "_trait")) %>%
        distinct %>%
        colMeans(., na.rm = TRUE),
      `Between SD` =
        data %>%
        select(paste0(selection, "_trait")) %>%
        distinct %>%
        summarise_all(sd, na.rm = TRUE) %>%
        unlist,
      `Within SD` = NA,
      `ICC(1)` =
        misty::multilevel.icc(
          data %>% select(selection),
          cluster = data[[id]],
          type = 1
        ),
      `ICC(2)` =
        misty::multilevel.icc(
          data %>% select(selection),
          cluster = data[[id]],
          type = 2
        )
    )
    
    
    for (i in 1:length(selection)) {
      dataRed <- 
        data %>%
        select(id, selection[i]) %>%
        na.omit
      descriptives$Within.SD[i] <-
        sqrt(sum(
          misty::cluster.scores(
            x = dataRed %>% select(selection[i]) ,
            cluster = dataRed %>% select(id),
            fun = "var",
            expand = FALSE
          ), 
          na.rm = TRUE
        ) / nrow(unique(data %>% select(id))))
    }
    
    descriptives <-
      descriptives %>%
      remove_rownames %>%
      column_to_rownames("variable") %>%
      round(2) %>%
      format(nsmall = 2)
    
    message("Upper Triangle: Between participants; Lower Triangle: Within participants")
    
    out <- 
      rbind(rMlComb, t(descriptives)) %>%
      as.data.frame
    rownames(out)[rownames(out) == "Grand.Mean"] <- "Grand Mean"
    rownames(out)[rownames(out) == "Between.SD"] <- "Between SD"
    rownames(out)[rownames(out) == "Within.SD"] <- "Within SD"
    rownames(out)[rownames(out) == "ICC.1."] <- "ICC(1)"
    rownames(out)[rownames(out) == "ICC.2."] <- "ICC(2)"
    
    return(out)
  }

MlTraitState <-
  function(data,
           id,
           selection) {
    
    # for testing
    # data <- dtWorkerSupp$workerInteractionType
    # data <- workerInteractionType
    # id <- "PID"
    # selection <- varSelection
    
    # calculate trait scores (participant means)
    data <- 
      data %>%
      #select(id, selection) %>%
      #mutate_if(is.factor, as.numeric) %>% # let's hope this doesn't cause any problems :-D
      mutate_at(selection, as.numeric) %>% # let's hope this works 
      mutate_at(selection, list(gm = ~mean(., na.rm=TRUE))) %>% # grand mean
      mutate_at(selection, list(gsd = ~sd(., na.rm=TRUE))) %>% # grand sd
      mutate_at(selection, list(gmc = ~.-mean(., na.rm=TRUE))) %>% # grand mean centered
      mutate_at(selection, list(gmz = ~(.-mean(., na.rm=TRUE))/sd(., na.rm=TRUE))) %>% # grand mean standardized
      group_by_at(vars(matches(id))) %>%
      mutate_at(selection, list(cm = ~mean(., na.rm=TRUE))) %>% # cluster mean
      mutate_at(selection, list(csd = ~sd(., na.rm=TRUE))) %>% # cluster sd
      ungroup 
    
    # center within cluster
    data <-
      map_dfc(selection,
              ~ data %>%
                transmute(
                  !!str_c(.x, '_cwc') :=
                    !!rlang::sym(.x)-!!rlang::sym(str_c(.x, "_cm"))
                )) %>%
      bind_cols(data, .)
    
    data <-
      map_dfc(selection,
              ~ data %>%
                transmute(
                  !!str_c(.x, '_zwc') :=
                    (!!rlang::sym(.x)-!!rlang::sym(str_c(.x, "_cm")))/!!rlang::sym(str_c(.x, "_csd"))
                )) %>%
      bind_cols(data, .)
    
    # cluster mean centered
    data <- 
      data %>%
      mutate_at(paste0(selection, "_cm"), list(c = ~.-mean(., na.rm=TRUE))) # grand mean
    
    data
  }


var_reduction = function(m0, m1) {
  library(tidyverse)
  lme4::VarCorr(m0) %>%
    as.data.frame %>%
    select(grp, var_m0 = vcov) %>%
    left_join(lme4::VarCorr(m1) %>%
                as.data.frame %>%
                select(grp, var_m1 = vcov)) %>%
    mutate(var_red = 1 - var_m1 / var_m0)
}

stdCoef.merMod <- function(object) {
  sdy <- sd(getME(object, "y"))
  sdx <- apply(getME(object, "X"), 2, sd)
  sc <- fixef(object) * sdx / sdy
  se.fixef <- coef(summary(object))[, "Std. Error"]
  se <- se.fixef * sdx / sdy
  return(data.frame(stdcoef = sc, stdse = se))
}

MlCoeffLatex <- function(lmeMdl = NULL, lmerCI = NULL, varName = NULL) {
  b <- coef(summary(lmeMdl))[varName,"Value"] %>% round(2) %>% format(nsmall=2)
  df <- lmeMdl$fixDF$X[varName] %>% as.numeric
  t <- coef(summary(lmeMdl))[varName,"t-value"] %>% round(2) %>% format(nsmall=2)
  p <- ifelse(coef(summary(lmeMdl))[varName,"p-value"]<.001, "< .001", paste0("= ", coef(summary(lmeMdl))[varName,"p-value"] %>% round(3) %>% format(nsmall=3)))
  
  if (is.null(lmerCI)) {
    paste0(
      "\\textit{b} = ", b, 
      ", \\textit{t}(", format(df, big.mark=","), ") = ", t,
      ", \\textit{p} ", p
    )
  } else {
    CIlwr <- lmerCI[varName, "2.5 %"] %>% round(2) %>% format(nsmall=2)
    CIupr <- lmerCI[varName, "97.5 %"] %>% round(2) %>% format(nsmall=2)
    
    paste0(
      "\\textit{b} = ", b, 
      ", \\textit{t}(", format(df, big.mark=","), ") = ", t,
      ", \\textit{p} ", p,
      ", \\textit{95\\%CI}[", CIlwr, ", ", CIupr, "]"
    )
  }
}

MlCoeffHtml <- function(lmeMdl = NULL, lmerCI = NULL, varName = NULL) {
  b <- coef(summary(lmeMdl))[varName,"Value"] %>% round(2) %>% format(nsmall=2)
  df <- lmeMdl$fixDF$X[varName] %>% as.numeric
  t <- coef(summary(lmeMdl))[varName,"t-value"] %>% round(2) %>% format(nsmall=2)
  p <- ifelse(coef(summary(lmeMdl))[varName,"p-value"]<.001, "< .001", paste0("= ", coef(summary(lmeMdl))[varName,"p-value"] %>% round(3) %>% format(nsmall=3)))
  
  if (is.null(lmerCI)) {
    paste0(
      "<i>b</i> = ", b, 
      ", <i>t</i>(", format(df, big.mark=","), ") = ", t,
      ", <i>p</i> ", p
    )
  } else {
    CIlwr <- lmerCI[varName, "2.5 %"] %>% round(2) %>% format(nsmall=2)
    CIupr <- lmerCI[varName, "97.5 %"] %>% round(2) %>% format(nsmall=2)
    
    paste0(
      "<i>b</i> = ", b, 
      ", <i>t</i>(", format(df, big.mark=","), ") = ", t,
      ", <i>p</i> ", p,
      ", <i>95%CI</i>[", CIlwr, ", ", CIupr, "]"
    )
  }
}

LmCoeffLatex <- function(lmMdl = NULL, varName = NULL, ci = TRUE) {
  # lmMdl <- lmAllAttFreqQualX
  # varName <- "SumContactNL_c"
  
  lmCI <-  confint(lmMdl)
  lmEta <- effectsize::eta_squared(lmMdl)
  
  b <- coef(summary(lmMdl))[varName,"Estimate"] %>% round(2) %>% format(nsmall=2)
  df <- lmMdl$df.residual %>% as.numeric
  t <- coef(summary(lmMdl))[varName,"t value"] %>% round(2) %>% format(nsmall=2)
  p <- ifelse(coef(summary(lmMdl))[varName,"Pr(>|t|)"]<.001, "< .001", paste0("= ", coef(summary(lmMdl))[varName,"Pr(>|t|)"] %>% round(3) %>% format(nsmall=3)))
  etaSqP <- lmEta$Eta2_partial[lmEta$Parameter == varName]
  
  if (ci == TRUE) {
    CIlwr <- lmCI[varName, "2.5 %"] %>% round(2) %>% format(nsmall=2)
    CIupr <- lmCI[varName, "97.5 %"] %>% round(2) %>% format(nsmall=2)
    
    paste0(
      "\\textit{b} = ", b, 
      ", \\textit{t}(", format(df, big.mark=","), ") = ", t,
      ", \\textit{p} ", p,
      ", \\textit{95\\%CI}[", CIlwr, ", ", CIupr, "]"
    )
  } else {
    paste0(
      "\\textit{b} = ", b, 
      ", \\textit{t}(", format(df, big.mark=","), ") = ", t,
      ", \\textit{p} ", p
    )
  }
}


linesep<-function(x,y=character()){
  if(!length(x))
    return(y)
  linesep(x[-length(x)], c(rep('',x[length(x)]-1),'\\addlinespace',y))  
}

select_best_lmer_model <- function(data, structure, dependent_var, prediction_form, return = "", control = list(opt = "nlmimb"), lmer_scale = FALSE) {
  # Set the warning option to -1 to suppress all warnings
  options(warn = -1)
  # lm data
  vars <- c(dependent_var, unique(unlist(strsplit(prediction_form, "\\W+"))), unique(unlist(strsplit(structure, "\\W+"))))
  data_lme <- data %>% 
    select(all_of(vars)) %>%
    na.omit
  
  # Define model formulas
  prediction_form_z <- gsub("(\\b[A-Za-z]*)C\\b", "\\1Z", prediction_form, perl = TRUE)
  dependent_var_z <- paste0(dependent_var, "Z")
  
  null_model_formula <- as.formula(paste(dependent_var, "~ 1"))
  random_model_formula <- as.formula(paste(dependent_var, "~", prediction_form))
  
  null_structure <- as.formula(paste("~ 1 |", structure))
  slope_structure <- as.formula(paste("~ 1 +", paste(unique(unlist(strsplit(prediction_form, "\\W+"))), collapse = " + "), "|", structure))
  
  null_lmer <- as.formula(paste0(dependent_var, " ~ 1 + (1 | ", structure, ")"))
  null_lmer_z <- as.formula(paste0(dependent_var_z, " ~ 1 + (1 | ", structure, ")"))
  intercept_lmer <- as.formula(paste0(dependent_var, " ~ ", prediction_form, " + (1 | ", structure, ")"))
  intercept_lmer_z <- as.formula(paste0(dependent_var_z, " ~ ", prediction_form_z, " + (1 | ", structure, ")"))
  slope_lmer <-
    as.formula(paste0(
      dependent_var, " ~ ", prediction_form," + (1 + ",
      paste(unique(unlist(strsplit(prediction_form, "\\W+"))), collapse = " + "),
      " | ", structure,")"
    ))
  slope_lmer_z <-
    as.formula(paste0(
      dependent_var_z, " ~ ", prediction_form_z," + (1 + ",
      paste(unique(unlist(strsplit(prediction_form_z, "\\W+"))), collapse = " + "),
      " | ", structure,")"
    ))
  
  # Fit the models
  
  # NULL MODEL
  null_model <- lme(
    fixed = null_model_formula,
    random = null_structure,
    data = data_lme,
    na.action = na.omit,
    control = control
  )
  lme_null_intervals <- intervals(null_model, which = "fixed")
  lme_null_ci <- lme_null_intervals$fixed %>% 
    as_tibble %>%
    mutate(coef = rownames(lme_null_intervals$fixed)) %>%
    select(coef, lower, upper)
  rm(lme_null_intervals)
  lme_coef_null <- 
    coef(summary(null_model)) %>% 
    as_tibble %>%
    mutate(coef = rownames(coef(summary(null_model)))) %>%
    left_join(lme_null_ci, by = "coef") %>%
    transmute(
      coef = coef,
      b = Value,
      se = Std.Error,
      df = DF,
      t = `t-value`,
      p = `p-value`,
      lwr = lower,
      upr = upper,
      coef_latex = paste0(
        "\\textit{b} = ", format(round(b, 2), nsmall = 2), 
        ", \\textit{t}(", format(df, big.mark=","), ") = ", format(round(t, 2), nsmall = 2),
        ", \\textit{p} ", ifelse(p < .001, "< .001", paste0("= ", format(round(p, 3), nsmall = 3))),
        ", \\textit{95\\%CI}[",  format(round(lwr, 2), nsmall = 2), ", ",  format(round(upr, 2), nsmall = 2), "]"
      ),
      coef_html = paste0(
        "<i>b</i> = ", format(round(b, 2), nsmall = 2), 
        ", <i>t</i>(", format(df, big.mark=","), ") = ", format(round(t, 2), nsmall = 2),
        ", <i>p</i> ", ifelse(p < .001, "< .001", paste0("= ", format(round(p, 3), nsmall = 3))),
        ", <i>95%CI</i>[",  format(round(lwr, 2), nsmall = 2), ", ",  format(round(upr, 2), nsmall = 2), "]"
      )
    )
  
  lmer_null <- lmer(null_lmer, data = data)
  lmer_null_ci <- confint(method = "Wald", lmer_null)
  lmer_null_z <- lmer(null_lmer_z, data = data)
  lmer_coef_null <- NULL
  
  
  # RANDOM INTERCEPTS MODEL
  random_intercept_model <- lme(
    fixed = random_model_formula,
    random = null_structure,
    data = data_lme,
    na.action = na.omit,
    control = control
  )
  lme_intercept_intervals <- intervals(random_intercept_model, which = "fixed")
  lme_intercept_ci <- lme_intercept_intervals$fixed %>% 
    as_tibble %>%
    mutate(coef = rownames(lme_intercept_intervals$fixed)) %>%
    select(coef, lower, upper)
  rm(lme_intercept_intervals)
  lme_coef_intercept <- 
    coef(summary(random_intercept_model)) %>% 
    as_tibble %>%
    mutate(coef = rownames(coef(summary(random_intercept_model)))) %>%
    left_join(lme_intercept_ci, by = "coef") %>%
    transmute(
      coef = coef,
      b = Value,
      se = Std.Error,
      df = DF,
      t = `t-value`,
      p = `p-value`,
      lwr = lower,
      upr = upper,
      coef_latex = paste0(
        "\\textit{b} = ", format(round(b, 2), nsmall = 2), 
        ", \\textit{t}(", format(df, big.mark=","), ") = ", format(round(t, 2), nsmall = 2),
        ", \\textit{p} ", ifelse(p < .001, "< .001", paste0("= ", format(round(p, 3), nsmall = 3))),
        ", \\textit{95\\%CI}[",  format(round(lwr, 2), nsmall = 2), ", ",  format(round(upr, 2), nsmall = 2), "]"
      ),
      coef_html = paste0(
        "<i>b</i> = ", format(round(b, 2), nsmall = 2), 
        ", <i>t</i>(", format(df, big.mark=","), ") = ", format(round(t, 2), nsmall = 2),
        ", <i>p</i> ", ifelse(p < .001, "< .001", paste0("= ", format(round(p, 3), nsmall = 3))),
        ", <i>95%CI</i>[",  format(round(lwr, 2), nsmall = 2), ", ",  format(round(upr, 2), nsmall = 2), "]"
      )
    )
  
  lmer_intercept <- lmer(intercept_lmer, data = data)
  lmer_intercept_ci <- confint(method = "Wald", lmer_intercept)
  lmer_intercept_z <- lmer(intercept_lmer_z, data = data)
  lmer_coef_intercept <- 
    summ(lmer_intercept, confint = TRUE, scale = lmer_scale)$coeftable %>%
    as_tibble %>%
    transmute(
      coef = rownames(summ(lmer_intercept)$coeftable),
      b = Est.,
      df = d.f.,
      t = `t val.`,
      p = p,
      lwr = `2.5%`,
      upr = `97.5%`,
      coef_latex = paste0(
        "\\textit{b} = ", format(round(b, 2), nsmall = 2), 
        ", \\textit{t}(", format(df, big.mark=","), ") = ", format(round(t, 2), nsmall = 2),
        ", \\textit{p} ", ifelse(p < .001, "< .001", paste0("= ", format(round(p, 3), nsmall = 3))),
        ", \\textit{95\\%CI}[",  format(round(lwr, 2), nsmall = 2), ", ",  format(round(upr, 2), nsmall = 2), "]"
      ),
      coef_html = paste0(
        "<i>b</i> = ", format(round(b, 2), nsmall = 2), 
        ", <i>t</i>(", format(df, big.mark=","), ") = ", format(round(t, 2), nsmall = 2),
        ", <i>p</i> ", ifelse(p < .001, "< .001", paste0("= ", format(round(p, 3), nsmall = 3))),
        ", <i>95%CI</i>[",  format(round(lwr, 2), nsmall = 2), ", ",  format(round(upr, 2), nsmall = 2), "]"
      )
    )
  
  
  # RANDOM SLOPES MODEL 
  random_slope_model <- lme(
    fixed = random_model_formula,
    random = slope_structure,
    data = data_lme,
    na.action = na.omit,
    control = control
  )
  lme_slope_intervals <- intervals(random_slope_model, which = "fixed")
  lme_slope_ci <- lme_slope_intervals$fixed %>% 
    as_tibble %>%
    mutate(coef = rownames(lme_slope_intervals$fixed)) %>%
    select(coef, lower, upper)
  rm(lme_slope_intervals)
  lme_coef_slope <- 
    coef(summary(random_slope_model)) %>% 
    as_tibble %>%
    mutate(coef = rownames(coef(summary(random_slope_model)))) %>%
    left_join(lme_slope_ci, by = "coef") %>%
    transmute(
      coef = coef,
      b = Value,
      se = Std.Error,
      df = DF,
      t = `t-value`,
      p = `p-value`,
      lwr = lower,
      upr = upper,
      coef_latex = paste0(
        "\\textit{b} = ", format(round(b, 2), nsmall = 2), 
        ", \\textit{t}(", format(df, big.mark=","), ") = ", format(round(t, 2), nsmall = 2),
        ", \\textit{p} ", ifelse(p < .001, "< .001", paste0("= ", format(round(p, 3), nsmall = 3))),
        ", \\textit{95\\%CI}[",  format(round(lwr, 2), nsmall = 2), ", ",  format(round(upr, 2), nsmall = 2), "]"
      ),
      coef_html = paste0(
        "<i>b</i> = ", format(round(b, 2), nsmall = 2), 
        ", <i>t</i>(", format(df, big.mark=","), ") = ", format(round(t, 2), nsmall = 2),
        ", <i>p</i> ", ifelse(p < .001, "< .001", paste0("= ", format(round(p, 3), nsmall = 3))),
        ", <i>95%CI</i>[",  format(round(lwr, 2), nsmall = 2), ", ",  format(round(upr, 2), nsmall = 2), "]"
      )
    )
  
  lmer_slope <- lmer(slope_lmer, data = data)
  lmer_slope_ci <- confint(method = "Wald", lmer_slope)
  lmer_slope_z <- lmer(slope_lmer_z, data = data)
  lmer_coef_slope <- 
    summ(lmer_slope, confint = TRUE, scale = lmer_scale)$coeftable %>%
    as_tibble %>%
    transmute(
      coef = rownames(summ(lmer_slope)$coeftable),
      b = Est.,
      df = d.f.,
      t = `t val.`,
      p = p,
      lwr = `2.5%`,
      upr = `97.5%`,
      coef_latex = paste0(
        "\\textit{b} = ", format(round(b, 2), nsmall = 2), 
        ", \\textit{t}(", format(df, big.mark=","), ") = ", format(round(t, 2), nsmall = 2),
        ", \\textit{p} ", ifelse(p < .001, "< .001", paste0("= ", format(round(p, 3), nsmall = 3))),
        ", \\textit{95\\%CI}[",  format(round(lwr, 2), nsmall = 2), ", ",  format(round(upr, 2), nsmall = 2), "]"
      ),
      coef_html = paste0(
        "<i>b</i> = ", format(round(b, 2), nsmall = 2), 
        ", <i>t</i>(", format(df, big.mark=","), ") = ", format(round(t, 2), nsmall = 2),
        ", <i>p</i> ", ifelse(p < .001, "< .001", paste0("= ", format(round(p, 3), nsmall = 3))),
        ", <i>95%CI</i>[",  format(round(lwr, 2), nsmall = 2), ", ",  format(round(upr, 2), nsmall = 2), "]"
      )
    )
  
  
  # Perform likelihood ratio tests
  anova_null_intercept <- anova(null_model, random_intercept_model)
  anova_intercept_slope <- anova(random_intercept_model, random_slope_model)
  anova_overall <- anova(null_model, random_intercept_model, random_slope_model)
  
  
  # If the random intercept model is significantly worse than the random slope model
  if ((anova_null_intercept$`p-value`[2] < 0.05 & anova_intercept_slope$`p-value`[2] < 0.05) | return == "slope") {
    return(
      list(
        choice = "random slopes model",
        anova = anova_overall,
        
        lme_model = random_slope_model,
        lme_ci = lme_slope_ci,
        lme_coef = lme_coef_slope,
        
        lmer_formula = gsub("\\s+", " ", paste(deparse(slope_lmer), collapse = "")),
        lmer_formula_z = gsub("\\s+", " ", paste(deparse(slope_lmer_z), collapse = "")),
        lmer_model = lmer_slope,
        lmer_ci = lmer_slope_ci,
        lmer_model_z = lmer_slope_z,
        lmer_coef = lmer_coef_slope
      )
    )
  } else if (anova_null_intercept$`p-value`[2] < 0.05 | return == "intercept") {
    return(
      list(
        choice = "random intercepts model",
        anova = anova_overall,
        
        lme_model = random_intercept_model,
        lme_ci = lme_intercept_ci,
        lme_coef = lme_coef_intercept,
        
        lmer_formula = gsub("\\s+", " ", paste(deparse(intercept_lmer), collapse = "")),
        lmer_formula_z = gsub("\\s+", " ", paste(deparse(intercept_lmer_z), collapse = "")),
        lmer_model = lmer_intercept,
        lmer_ci = lmer_intercept_ci,
        lmer_model_z = lmer_intercept_z,
        lmer_coef = lmer_coef_intercept
      )
    )
  } else {
    return(
      list(
        choice = "null model",
        anova = anova_overall,
        
        lme_model = null_model,
        lme_ci = lme_null_ci,
        lme_coef = lme_coef_null,
        
        lmer_formula = gsub("\\s+", " ", paste(deparse(null_lmer), collapse = "")),
        lmer_formula_z = gsub("\\s+", " ", paste(deparse(null_lmer_z), collapse = "")),
        lmer_model = lmer_null,
        lmer_ci = lmer_null_ci,
        lmer_model_z = lmer_null_z,
        lmer_coef = lmer_coef_null
      )
    )
  }
}

get_latex_coef <- function(data, coef_name) {
  latex_out <- data[data$coef == coef_name, "coef_latex"]
  return(latex_out)
}
