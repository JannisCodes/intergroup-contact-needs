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
      `Between-person SD` =
        data %>%
        select(paste0(selection, "_trait")) %>%
        distinct %>%
        summarise_all(sd, na.rm = TRUE) %>%
        unlist,
      `Within-person SD` = NA,
      `ICC(1)` =
        misty::multilevel.icc(
          data %>% select(selection),
          cluster = data$PID,
          type = 1
        ),
      `ICC(2)` =
        misty::multilevel.icc(
          data %>% select(selection),
          cluster = data$PID,
          type = 2
        )
    )
    
    for (i in 1:length(selection)) {
      dataRed <- 
        data %>%
        select(id, selection[i]) %>%
        na.omit
      descriptives$Within.person.SD[i] <-
        sqrt(sum(
          misty::cluster.scores(
            x = dataRed %>% select(selection[i]) ,
            cluster = dataRed %>% select(id),
            fun = "var",
            expand = FALSE
          )
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
    rownames(out)[rownames(out) == "Between.person.SD"] <- "Between-person SD"
    rownames(out)[rownames(out) == "Within.person.SD"] <- "Within-person SD"
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
    # id <- "PID"
    # selection <- varSelection
    
    # calculate trait scores (participant means)
    data <- 
      data %>%
      #select(id, selection) %>%
      mutate_at(selection, list(gm = ~mean(., na.rm=TRUE))) %>% # grand mean
      mutate_at(selection, list(gmc = ~.-mean(., na.rm=TRUE))) %>% # grand mean centered
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


var_reduction = function(m0, m1){
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
      ", t(", df, ") = ", t,
      ", \\textit{p} = ", p
    )
  } else {
    CIlwr <- lmerCI[varName, "2.5 %"] %>% round(2) %>% format(nsmall=2)
    CIupr <- lmerCI[varName, "97.5 %"] %>% round(2) %>% format(nsmall=2)
    
    paste0(
      "\\textit{b} = ", b, 
      ", t(", df, ") = ", t,
      ", \\textit{p} = ", p,
      ", \\textit{95%CI[", CIlwr, ", ", CIupr, "]"
    )
  }
}
