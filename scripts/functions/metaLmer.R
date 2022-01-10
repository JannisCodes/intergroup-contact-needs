metaLmerOut <- function(lmerDataTbl, type = "FE", name, title, ...) {
  # lmerDataTbl <- theoryMetaTbl
  # name <- "Theory"
  # title <- "Need Fullfillment — Mixed Effects Regression"
  
  # extract parameter names
  var <-
    as.character(unique(lmerDataTbl$coef))[!as.character(unique(lmerDataTbl$coef)) %in% c("(Intercept)")]
  names <-
    stri_replace_all_fixed(
      var,
      pattern = c(
        "SumContactNL", 
        "AvQuality", 
        "SumContactNL.AvQuality",
        "NonOutgroupInteraction",
        "OutgroupInteraction",
        "CoreNeedZ",
        "QualityZ",
        "_zwc"
      ),
      replacement = c(
        "Number of Outgroup Contacts",
        "Average Interaction Quality",
        "Contact X Quality",
        "Non-Outgroup Interaction",
        "Outgroup Interaction",
        "Core Need Fulfillment",
        "Interaction Quality",
        ""
      ),
      vectorize_all = FALSE
    )
  studies <- as.character(unique(lmerDataTbl$sample))
  
  ### Parametric Meta-Analysis
  meta.parametric <- list()
  
  sum.names <- c("effect", 
                 "name",
                 paste(c("beta", "CI"),rep(studies, each = 2),sep="."),
                 "meta.beta",
                 "meta.CI",
                 "meta.pval",
                 "meta.star")
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
      var.names = c(paste(var[i], ".yi", sep = ""), paste(var[i], ".vi", sep =
                                                            ""))
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
    meta.parametric.sum$meta.beta[i] <- meta.parametric[[var[i]]]$b
    meta.parametric.sum$meta.CI[i] <- paste(
      "[",
      round(meta.parametric[[var[i]]]$ci.lb, 2),
      ", ",
      round(meta.parametric[[var[i]]]$ci.ub, 2),
      "]",
      sep = ""
    )
    meta.parametric.sum$meta.pval[i] <- meta.parametric[[var[i]]]$pval
    meta.parametric.sum$meta.star[i] <-
      ifelse(
        meta.parametric.sum$meta.pval[i] < 0.001,
        "***",
        ifelse(
          meta.parametric.sum$meta.pval[i] < 0.01,
          "**",
          ifelse(meta.parametric.sum$meta.pval[i] < 0.05, "*", "")
        )
      )
    
    betaCI <- cbind(data.frame(effect = var[i]),
          data.frame(t(setNames(
            as.numeric(meta.parametric[[var[i]]]$yi),
            paste("beta", attributes(meta.parametric[[var[i]]]$yi)$slab, sep = ".")
          ))),
          data.frame(t(setNames(
            paste("[",
                  round(as.numeric(
                    as.numeric(meta.parametric[[var[i]]]$yi) -
                      1.96 * sqrt(as.numeric(meta.parametric[[var[i]]]$vi))
                  ), 2),
                  ", ",
                  round(as.numeric(
                    as.numeric(meta.parametric[[var[i]]]$yi) + 1.96 * sqrt(as.numeric(meta.parametric[[var[i]]]$vi))
                  ), 2),
                  "]",
                  sep = ""),
            paste("CI", attributes(meta.parametric[[var[i]]]$yi)$slab, sep = ".")
          ))))
    
    library(rqdatatable)
    meta.parametric.sum <- natural_join(meta.parametric.sum,
                              betaCI,
                              by = "effect",
                              jointype = "FULL") %>%
      select(names(meta.parametric.sum))
  }
  #meta.parametric
  
  eff <- escalc(
    data = lmerDataTbl,
    slab = coef,
    measure = "SPCOR",
    ti = tval,
    ni = n,
    mi = m,
    r2i = Rsq
  )
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
      #xlim = c(-2, 2),
      cex = 0.75,
      ylim = c(-1, max(rows)+4),
      rows = rows,
      slab = rep(c("Study 1", "Study 2", "Study 3"), length(var)),
      mlab = "",
      main = paste0("Meta Analysis: Forest Plot [Parametric] \n(", title, ")")
    ) 
    par(font = 2)
    text(pltForest$xlim[1], c(11, 5), cex = 0.75, pos = 4, names)
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
    meta.bootstrapped.sum$meta.beta[i] <- meta.bootstrapped[[var[i]]]$b
    meta.bootstrapped.sum$meta.CI[i] <- paste(
      "[",
      round(meta.bootstrapped[[var[i]]]$ci.lb, 2),
      ", ",
      round(meta.bootstrapped[[var[i]]]$ci.ub, 2),
      "]",
      sep = ""
    )
    meta.bootstrapped.sum$meta.pval[i] <- meta.bootstrapped[[var[i]]]$pval
    meta.bootstrapped.sum$meta.star[i] <-
      ifelse(
        meta.bootstrapped.sum$meta.pval[i] < 0.001,
        "***",
        ifelse(
          meta.bootstrapped.sum$meta.pval[i] < 0.01,
          "**",
          ifelse(meta.bootstrapped.sum$meta.pval[i] < 0.05, "*", "")
        )
      )
    
    betaCI <- cbind(data.frame(effect = var[i]),
                    data.frame(t(setNames(
                      as.numeric(meta.bootstrapped[[var[i]]]$yi),
                      paste("beta", attributes(meta.bootstrapped[[var[i]]]$yi)$slab, sep = ".")
                    ))),
                    data.frame(t(setNames(
                      paste("[",
                            round(as.numeric(
                              as.numeric(meta.bootstrapped[[var[i]]]$yi) -
                                1.96 * sqrt(as.numeric(meta.bootstrapped[[var[i]]]$vi))
                            ), 2),
                            ", ",
                            round(as.numeric(
                              as.numeric(meta.bootstrapped[[var[i]]]$yi) + 1.96 * sqrt(as.numeric(meta.bootstrapped[[var[i]]]$vi))
                            ), 2),
                            "]",
                            sep = ""),
                      paste("CI", attributes(meta.bootstrapped[[var[i]]]$yi)$slab, sep = ".")
                    ))))
    
    meta.bootstrapped.sum <- natural_join(meta.bootstrapped.sum,
                                          betaCI,
                                          by = "effect",
                                          jointype = "FULL") %>%
      select(names(meta.bootstrapped.sum))
  }
  #meta.bootstrapped
  
  eff.boot <-
    escalc(
      data = lmerDataTbl %>% filter(coef != "(Intercept)"),
      slab = coef,
      measure = "SPCOR",
      yi = B,
      sei = se,
      ni = n,
      mi = m
    )
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
      main = paste0("Meta Analysis: Forest Plot [Parametric] \n(", title, ")")
    ) 
    metafor::forest(
      meta.boot,
      #xlim = c(-2, 2),
      cex = 0.75,
      ylim = c(-1, max(rows)+4),
      rows = rows,
      slab = rep(c("Study 1", "Study 2", "Study 3"), length(var)),
      mlab = "",
      main = paste0("Meta Analysis: Forest Plot [Parametric] \n(", title, ")")
    ) 
    par(font = 2)
    text(pltForest$xlim[1], c(11, 5), cex = 0.75, pos = 4, names)
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
  
  list(tbl.parametric = meta.parametric.sum, tbl.bootstrapped = meta.bootstrapped.sum, meta.parametric = meta.parametric, meta.boot = meta.bootstrapped)
}
