strLabFix <- function(varIn, ...) {
  stri_replace_all_fixed(
    varIn,
    pattern = c(
      "SumContactNL", 
      "AvQuality", 
      "SumContactNL.AvQuality",
      "NonOutgroupInteraction",
      "OutgroupInteraction",
      "CoreNeedZ",
      "CoreNeed",
      "QualityZ",
      "_zwc",
      "."
    ),
    replacement = c(
      "Number of Outgroup Contacts",
      "Average Interaction Quality",
      "Contact X Quality",
      "Non-Outgroup Interaction",
      "Outgroup Interaction",
      "Core Need Fulfillment",
      "Core Need Fulfillment",
      "Interaction Quality",
      "",
      " X "
    ),
    vectorize_all = FALSE
  )
}

strStudyFix <- function(studyIn, ...) {
  stri_replace_all_fixed(
    studyIn,
    pattern = c(
      "Worker", 
      "worker", 
      "Student",
      "student",
      "Medical",
      "medical"
    ),
    replacement = c(
      "Study 1",
      "Study 1",
      "Study 2",
      "Study 2",
      "Study 3",
      "Study 3"
    ),
    vectorize_all = FALSE
  )
}

metaLmerOut <- function(lmerDataTbl, type = "FE", name, title, ...) {
  # lmerDataTbl <- GeneralLmerMetaTbl
  # name <- "Theory"
  # title <- "Need Fullfillment — Mixed Effects Regression"
  
  # extract parameter names
  var <- as.character(unique(lmerDataTbl$coef))[!as.character(unique(lmerDataTbl$coef)) %in% c("(Intercept)")]
  names <- strLabFix(var)
  studies <- as.character(unique(lmerDataTbl$sample))
  
  
  ### Parametric Meta-Analysis
  meta.parametric <- list()
  
  sum.names <- c("effect", 
                 "name",
                 paste(c("beta", "lwr", "upr", "CI"),rep(studies, each = 4),sep="."),
                 "beta.meta",
                 "lwr.meta",
                 "upr.meta",
                 "CI.meta",
                 "pval.meta",
                 "star.meta")
  meta.parametric.sum <-
    data.frame(matrix(
      ncol = length(sum.names),
      nrow = length(var),
      dimnames = list(NULL, sum.names)
    )) %>%
    mutate(
      effect = var,
      name = names
    )
  
  # i = 1
  for (i in 1:length(var)) {
    md <- escalc(
      data = lmerDataTbl[lmerDataTbl$coef == var[i], ],
      measure = "SPCOR",
      ti = tval,
      ni = n,
      mi = m,
      r2i = Rsq,
      var.names = c(paste(var[i], ".yi", sep = ""), paste(var[i], ".vi", sep = ""))
    )
    
    meta.parametric[[var[i]]] <-
      rma(
        yi = get(paste(var[i], ".yi", sep = "")),
        vi = get(paste(var[i], ".vi", sep =
                         "")),
        data = md,
        method = type,
        slab = sample
      )
    
    meta.parametric.sum$effect[i] <- var[i]
    meta.parametric.sum$beta.meta[i] <- meta.parametric[[var[i]]]$b
    meta.parametric.sum$lwr.meta[i] <- meta.parametric[[var[i]]]$ci.lb
    meta.parametric.sum$upr.meta[i] <- meta.parametric[[var[i]]]$ci.ub
    meta.parametric.sum$CI.meta[i] <- paste(
      "[",
      format(round(meta.parametric[[var[i]]]$ci.lb, 2), nsmall = 2),
      ", ",
      format(round(meta.parametric[[var[i]]]$ci.ub, 2), nsmall = 2),
      "]",
      sep = ""
    )
    meta.parametric.sum$pval.meta[i] <- meta.parametric[[var[i]]]$pval
    meta.parametric.sum$star.meta[i] <-
      ifelse(
        meta.parametric.sum$pval.meta[i] < 0.001,
        "***",
        ifelse(
          meta.parametric.sum$pval.meta[i] < 0.01,
          "**",
          ifelse(meta.parametric.sum$pval.meta[i] < 0.05, "*", "")
        )
      )
    
    betaCI <- cbind(data.frame(effect = var[i]),
          data.frame(t(setNames(
            as.numeric(meta.parametric[[var[i]]]$yi),
            paste("beta", attributes(meta.parametric[[var[i]]]$yi)$slab, sep = ".")
          ))),
          data.frame(t(setNames(
            as.numeric(as.numeric(meta.parametric[[var[i]]]$yi) -
                         1.96 * sqrt(as.numeric(
                           meta.parametric[[var[i]]]$vi
                         ))),
            paste("lwr", attributes(meta.parametric[[var[i]]]$yi)$slab, sep = ".")
          ))), 
          data.frame(t(setNames(
            as.numeric(as.numeric(meta.parametric[[var[i]]]$yi) +
                         1.96 * sqrt(as.numeric(
                           meta.parametric[[var[i]]]$vi
                         ))),
            paste("upr", attributes(meta.parametric[[var[i]]]$yi)$slab, sep = ".")
          ))), 
          data.frame(t(setNames(
            paste("[",
                  format(round(as.numeric(
                    as.numeric(meta.parametric[[var[i]]]$yi) -
                      1.96 * sqrt(as.numeric(meta.parametric[[var[i]]]$vi))
                  ), 2), nsmall = 2),
                  ", ",
                  format(round(as.numeric(
                    as.numeric(meta.parametric[[var[i]]]$yi) + 1.96 * sqrt(as.numeric(meta.parametric[[var[i]]]$vi))
                  ), 2), nsmall = 2),
                  "]",
                  sep = ""),
            paste("CI", attributes(meta.parametric[[var[i]]]$yi)$slab, sep = ".")
          ))))
    
    for(r in 2:ncol(betaCI)) {
      meta.parametric.sum[meta.parametric.sum$effect==betaCI$effect, names(betaCI[r])] <- betaCI[r]
    }
  }
  #meta.parametric
  #meta.parametric.sum
  
  eff <- escalc(
    data = lmerDataTbl,
    slab = coef,
    measure = "SPCOR",
    ti = tval,
    ni = n,
    mi = m,
    r2i = Rsq
  ) %>%
    mutate(effLwr = yi - 1.96 * sqrt(vi),
           effUpr = yi + 1.96 * sqrt(vi),
           effCI = paste("[", format(round(effLwr, 2), nsmall = 2), ", ", format(round(effUpr, 2), nsmall = 2), "]", sep = ""),
           name = strLabFix(coef),
           study = strStudyFix(sample))
    
  meta <- rma(
    yi = yi,
    vi = vi,
    data = eff,
    method = type,
    slab = paste(coef, sample, sep = " - ")
  )
  
  k <- length(var) * length(studies)
  rows <- c()
  for (j in seq(0, by = length(studies)+3, length.out = length(var))) {
    rowsTmp <- seq(2+j, length.out = length(studies))
    rows <- c(rows, rowsTmp)
  }
  rows <- rev(rows)
  png(file = paste0("Figures/forestParametric", type, name, ".png"), height = 200 + 40 * k ^ .85)
    # save details for text placement
    pltForest <- metafor::forest(
      meta,
      #xlim = c(-1, 1),
      cex = 0.75,
      ylim = c(-1, max(rows)+4),
      rows = rows,
      slab = rep(c("Study 1", "Study 2", "Study 3"), length(var)),
      mlab = "",
      main = paste0("Meta Analysis: Forest Plot [Parametric] \n(", title, ")"),
      annotate = TRUE, 
      addfit = FALSE, 
      addpred = FALSE,
      showweights = FALSE, 
      header = FALSE
    ) 
    par(font = 2)
    text(pltForest$xlim[1], rev(seq(length(studies)+2, by = length(studies)+3, length.out = length(var))), cex = 0.75, pos = 4, names)
    m = 1
    for (i in length(var):1) {
      addpoly(
        meta.parametric[[i]],
        row = m,
        cex = 0.75,
        mlab = paste("►", mapvalues(
          type,
          from = c(
            "FE",
            "REML"
          ),
          to = c(
            "Fixed Effect",
            "Random Effect"
          )))
      )
      m <- m + 6
    }
    text(pltForest$xlim[1],
         max(rows)+3,
         "Predictors by study",
         cex = .8,
         pos = 4)
    text(pltForest$xlim[2],
         max(rows)+3,
         "Semi-Partial Correlation [95% CI]",
         cex = .8,
         pos = 2)
  dev.off()
  
  ### Bootstrapped Meta-Analysis
  var.boot <- var
  names.boot <- names
  meta.bootstrapped <- list()
  meta.bootstrapped.sum <-
    data.frame(matrix(
      ncol = length(sum.names),
      nrow = length(var),
      dimnames = list(NULL, sum.names)
    )) %>%
    mutate(
      effect = var.boot,
      name = names.boot
    )
  
  # i = 1
  for (i in 1:length(var.boot)) {
    md <- escalc(
      data = lmerDataTbl[lmerDataTbl$coef == var.boot[i], ],
      measure = "SPCOR",
      yi = B,
      sei = se,
      ni = n,
      mi = m,
      var.names = c(
        paste(var.boot[i], ".yi", sep = ""),
        paste(var.boot[i], ".vi", sep = "")
      )
    )
    
    meta.bootstrapped[[var.boot[i]]] <-
      rma(
        yi = get(paste(var.boot[i], ".yi", sep = "")),
        vi = get(paste(var.boot[i], ".vi", sep =
                         "")),
        data = md,
        method = "FE",
        slab = sample
      )
    
    meta.bootstrapped.sum$effect[i] <- var[i]
    meta.bootstrapped.sum$beta.meta[i] <- meta.bootstrapped[[var[i]]]$b
    meta.bootstrapped.sum$lwr.meta[i] <- meta.bootstrapped[[var[i]]]$ci.lb
    meta.bootstrapped.sum$upr.meta[i] <- meta.bootstrapped[[var[i]]]$ci.ub
    meta.bootstrapped.sum$CI.meta[i] <- paste(
      "[",
      format(round(meta.bootstrapped[[var[i]]]$ci.lb, 2), nsmall = 2),
      ", ",
      format(round(meta.bootstrapped[[var[i]]]$ci.ub, 2), nsmall = 2),
      "]",
      sep = ""
    )
    meta.bootstrapped.sum$pval.meta[i] <- meta.bootstrapped[[var[i]]]$pval
    meta.bootstrapped.sum$star.meta[i] <-
      ifelse(
        meta.bootstrapped.sum$pval.meta[i] < 0.001,
        "***",
        ifelse(
          meta.bootstrapped.sum$pval.meta[i] < 0.01,
          "**",
          ifelse(meta.bootstrapped.sum$pval.meta[i] < 0.05, "*", "")
        )
      )
    
    betaCI <- cbind(data.frame(effect = var[i]),
                    data.frame(t(setNames(
                      as.numeric(meta.bootstrapped[[var[i]]]$yi),
                      paste("beta", attributes(meta.bootstrapped[[var[i]]]$yi)$slab, sep = ".")
                    ))),
                    data.frame(t(setNames(
                      as.numeric(as.numeric(meta.bootstrapped[[var[i]]]$yi) -
                                   1.96 * sqrt(as.numeric(
                                     meta.bootstrapped[[var[i]]]$vi
                                   ))),
                      paste("lwr", attributes(meta.bootstrapped[[var[i]]]$yi)$slab, sep = ".")
                    ))), 
                    data.frame(t(setNames(
                      as.numeric(as.numeric(meta.bootstrapped[[var[i]]]$yi) +
                                   1.96 * sqrt(as.numeric(
                                     meta.bootstrapped[[var[i]]]$vi
                                   ))),
                      paste("upr", attributes(meta.bootstrapped[[var[i]]]$yi)$slab, sep = ".")
                    ))), 
                    data.frame(t(setNames(
                      paste("[",
                            format(round(as.numeric(
                              as.numeric(meta.bootstrapped[[var[i]]]$yi) -
                                1.96 * sqrt(as.numeric(meta.bootstrapped[[var[i]]]$vi))
                            ), 2), nsmall = 2),
                            ", ",
                            format(round(as.numeric(
                              as.numeric(meta.bootstrapped[[var[i]]]$yi) + 1.96 * sqrt(as.numeric(meta.bootstrapped[[var[i]]]$vi))
                            ), 2), nsmall = 2),
                            "]",
                            sep = ""),
                      paste("CI", attributes(meta.bootstrapped[[var[i]]]$yi)$slab, sep = ".")
                    ))))
    for(r in 2:ncol(betaCI)) {
      meta.bootstrapped.sum[meta.bootstrapped.sum$effect==betaCI$effect, names(betaCI[r])] <- betaCI[r]
    }
  }
  #meta.bootstrapped
  #meta.bootstrapped.sum
  
  eff.boot <-
    escalc(
      data = lmerDataTbl %>% filter(coef != "(Intercept)"),
      slab = coef,
      measure = "SPCOR",
      yi = B,
      sei = se,
      ni = n,
      mi = m
    ) %>%
    mutate(effLwr = yi - 1.96 * sqrt(vi),
           effUpr = yi + 1.96 * sqrt(vi),
           effCI = paste("[", format(round(effLwr, 2), nsmall = 2), ", ", format(round(effUpr, 2), nsmall = 2), "]", sep = ""),
           name = strLabFix(coef),
           study = strStudyFix(sample))
  
  meta.boot <- rma(
    yi = yi,
    vi = vi,
    data = eff.boot,
    method = type,
    slab = paste(coef, sample, sep = " - ")
  )
  
  png(file = paste0("Figures/forestBootstrapped", type, name, ".png"), height = 200 + 40 * k ^ .85)
    # save details for text placement
    pltForest <- metafor::forest(
      meta.boot,
      #xlim = c(-2, 2),
      cex = 0.75,
      ylim = c(-1, max(rows)+4),
      rows = rows,
      slab = rep(c("Study 1", "Study 2", "Study 3"), length(var)),
      mlab = "",
      main = paste0("Meta Analysis: Forest Plot [Parametric] \n(", title, ")"),
      annotate = TRUE, 
      addfit = FALSE, 
      addpred = FALSE,
      showweights = FALSE, 
      header = FALSE
    ) 
    par(font = 2)
    text(pltForest$xlim[1], rev(seq(length(studies)+2, by = length(studies)+3, length.out = length(var))), cex = 0.75, pos = 4, names)
    m = 1
    for (i in length(var):1) {
      addpoly(
        meta.parametric[[i]],
        row = m,
        cex = 0.75,
        mlab = paste("►", mapvalues(
          type,
          from = c(
            "FE",
            "REML"
          ),
          to = c(
            "Fixed Effect",
            "Random Effect"
          )))
      )
      m <- m + 6
    }
    text(pltForest$xlim[1],
         max(rows)+3,
         "Predictors by study",
         cex = .8,
         pos = 4)
    text(pltForest$xlim[2],
         max(rows)+3,
         "Semi-Partial Correlation [95% CI]",
         cex = .8,
         pos = 2)
  dev.off()
  
  tbl.parametric <- meta.parametric.sum %>%
    select(-starts_with("lwr"),
           -starts_with("upr"))
  
  tbl.bootstrapped <- meta.bootstrapped.sum %>%
    select(-starts_with("lwr"),
           -starts_with("upr"))
  
  list(tbl.parametric = tbl.parametric, meta.parametric.sum = meta.parametric.sum, eff = eff, meta.parametric = meta.parametric, 
       tbl.bootstrapped = tbl.bootstrapped, meta.bootstrapped.sum = meta.bootstrapped.sum, eff.boot = eff.boot, meta.boot = meta.boot)
}










metaSubForest <- function(effects, title = "", addAbove = 0, filename = NULL, width = 600, height = 800,...) {
  metaGeneral <- rma(
    yi = yi,
    vi = vi,
    data = effects,
    method = "FE",
    slab = paste(analysis, coef, sample, sep = " - ")
  )
  
  
  m <- 0
  studyIndex <- list()
  varIndex <- list()
  metaIndex <- list()
  analysisIndex <- list()
  for (i in length(unique(effects$analysis)):1) {
    studyIndex[[i]] <- as.integer(sapply(
      seq(
        from = m,
        by = 6,
        length.out = length(unique(effects$coef[effects$analysis == unique(effects$analysis)[i]]))
      ),
      seq,
      length.out = length(unique(effects$sample[effects$analysis == unique(effects$analysis)[i]]))
    ))
    
    varIndex[[i]] <- seq(
      from = m + length(unique(effects$sample[effects$analysis == unique(effects$analysis)[i]])),
      by = 6,
      length.out = length(unique(effects$coef[effects$analysis == unique(effects$analysis)[i]]))
    )
    
    metaIndex[[i]] <- seq(
      from = m - 1,
      by = 6,
      length.out = length(unique(effects$coef[effects$analysis == unique(effects$analysis)[i]]))
    )
    
    analysisIndex[[i]] <- max(varIndex[[i]]) + 1
    
    m <- max(studyIndex[[i]]) + 5
  }
  
  
  subgroupMeta <- list()
  for (i in length(unique(effects$metaId)):1){
    subgroupMeta[[i]] <- rma(yi, 
                             vi, 
                             data = effects %>% filter(metaId == unique(effects$metaId)[i]), 
                             method = "FE",
                             slab = paste(analysis, coef, sep = " - "))
  }
  
  png(file = paste0("Figures/forestParametric", filename, ".png"), width = width, height = height)
  forestGeneral <- forest(
    metaGeneral,
    cex = 0.75,
    ylim = c(-1, max(unlist(studyIndex)) + addAbove),
    # alim = c(-1, 1), # only for General Pt.1
    # xlim = c(-2, 2), # only for General Pt.1
    # xlab = NULL, # only for General Pt.1
    # alim = c(-0.2, 0.2), # only for General Pt.2
    # xlim = c(-0.4, 0.4), # only for General Pt.2
    #order = effects$analysis,
    slab = effects$study,
    rows = rev(unlist(rev(studyIndex))),
    main = title,
    annotate = TRUE,
    addfit = FALSE,
    addpred = FALSE,
    showweights = FALSE,
    header = FALSE
  )
  text(forestGeneral$xlim[1],
       forestGeneral$ylim[2]-1,
       "Predictors by study",
       cex = .8,
       font = 2,
       pos = 4)
  text(
    forestGeneral$xlim[2],
    forestGeneral$ylim[2]-1,
    "Semi-Partial Correlation [95% CI]",
    cex = .8,
    font = 2,
    pos = 2
  )
  text(
    forestGeneral$xlim[1], 
    rev(unlist(rev(varIndex))), 
    #unique(paste(effects$analysis, effects$coef, sep = " - ")), # for order checking
    # stri_replace_all_regex(unique(paste(effects$analysis, effects$name, sep = " - ")),
    #                        pattern = paste0(unique(effects$analysis), " - "),
    #                        replacement=rep("", length(unique(effects$analysis))),
    #                        vectorize=FALSE),
    gsub(".* - " , "", unique(paste(effects$analysis, effects$name, sep = " - "))),
    cex = .8,
    font = 4,
    pos = 4
  )
  text(
    forestGeneral$xlim[1], 
    rev(unlist(rev(analysisIndex))), 
    unique(effects$analysis),
    cex = .8,
    font = 2,
    pos = 4
  )
  par(font = 3)
  for (i in 1:length(unique(effects$metaId))) {
    addpoly(
      subgroupMeta[[i]], 
      row = rev(unlist(rev(metaIndex)))[[i]],
      mlab = "► Fixed Effect",
      cex = 0.75
    )
  }
  dev.off()
  
}

