kappa.full <- function(x, name, ...){
  xtab <- table(x)
  diagonal.counts <- diag(xtab)
  N <- sum(xtab)
  row.marginal.props <- rowSums(xtab)/N
  col.marginal.props <- colSums(xtab)/N
  # Compute kappa (k)
  Po <- sum(diagonal.counts)/N
  Pe <- sum(row.marginal.props*col.marginal.props)
  k <- (Po - Pe)/(1 - Pe)
  se <- sqrt((Po*(1-Po))/(N*(1-Pe)^2)) # according to PMID: 23092060
  lwr <- k-1.96*se
  upr <- k+1.96*se
  
  data.frame(name,Po,Pe,k,se,lwr,upr)
}

# # Minimally reproducible dataset
# coders <- data.table(Affect1 = rbinom(250, 1,.62), Affect2 = rbinom(250, 1,.58))
# name <- "affect"
# kappa.full(coders, name)

# --------------------------------------------------------------------------------

kappa.full.multiple <- function(x, names, ...){
  kComb <- list()
  for (i in 1:length(names)) {
    kComb[[i]] <- kappa.full(x %>% select(contains(names[i])), names[i])
  }
  do.call(rbind.data.frame, kComb)  
}

# # Minimally reproducible dataset
# coders <- data.table(Affect1 = rbinom(250, 1,.62), Affect2 = rbinom(250, 1,.58),
#                      Behavior1 = rbinom(250, 1,.32), Behavior2 = rbinom(250, 1,.40),
#                      Cognition1 = rbinom(250, 1,.52), Cognition2 = rbinom(250, 1,.55),
#                      Desire1 = rbinom(250, 1,.78), Desire2 = rbinom(250, 1,.82))
# names <- c("affect", "behavior", "cognition", "desire")
# kappa.full.multiple(coders, names)

# --------------------------------------------------------------------------------

kappa.pooled <- function(ktab, ...){
  (mean(ktab$Po, na.rm = TRUE) - mean(ktab$Pe, na.rm = TRUE))/(1 - mean(ktab$Pe, na.rm = TRUE))
}

kappa.pooled2 <- function(ktab, ...){
  k <- kappa.pooled(ktab)
  require(boot)
  set.seed(42)
  boot <- boot(ktab, function(x, idx) {kappa.pooled(x[idx,])}, R = 10000)
  se <- sd(boot$t)
  ci <- boot.ci(boot.out = boot, type = "basic")
  lwr <- ci$basic[4]
  upr <- ci$basic[5]
  
  data.frame(k.pooled=k, se, lwr, upr)
}

# # Minimally reproducible dataset
# coders <- data.table(Affect1 = rbinom(250, 1,.62), Affect2 = rbinom(250, 1,.58),
#                      Behavior1 = rbinom(250, 1,.32), Behavior2 = rbinom(250, 1,.40),
#                      Cognition1 = rbinom(250, 1,.52), Cognition2 = rbinom(250, 1,.55),
#                      Desire1 = rbinom(250, 1,.78), Desire2 = rbinom(250, 1,.82))
# names <- c("affect", "behavior", "cognition", "desire")
# dt.kappa <- kappa.full.multiple(coders, names)
# kappa.pooled2(dt.kappa, names)
