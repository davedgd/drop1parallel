library(doParallel)
library(MASS)
library(gtools)
library(lme4)

drop1.parallel <- function (theModel, test = "Chisq") { # only supports Chisq curently for all models (argument not used)
  
  if (class(theModel) != "glmerMod")
    stop("You must supply a glmer model...")

  originalData <- theModel@call[['data']]
  
  if (length(intersect(class(theModel), c("glmerMod", "lmerMod", "lmerModLmerTest"))) > 0) {
    potentialDrops <- drop.scope(lme4:::getFixedFormula(formula(theModel)), formula("~ 1"))
    cat(paste("Dropping ", length(potentialDrops), " terms on ", getDoParWorkers(), " cores...\n\n", sep = ""))
    res <- foreach(theDrop = potentialDrops, .packages = attr(class(theModel), "package")) %dopar% update(theModel, as.formula(paste(".~.-", theDrop)), data = theModel@frame) # data = dat
  } else {
    potentialDrops <- drop.scope(formula(theModel), formula("~ 1"))
    cat(paste("Dropping ", length(potentialDrops), " terms on ", getDoParWorkers(), " cores...\n\n", sep = ""))
    res <- foreach(theDrop = potentialDrops, .packages = attr(class(theModel), "package")) %dopar% update(theModel, as.formula(paste(".~.-", theDrop)))
  }
  
  for (each in 1:length(res))
    res[[each]]@call[['data']] <- originalData
  #lapply(res, function(x) x@call[['data']])
  #anova(theModel, res[[1]])
  
  cat("Single term deletions\n\nModel:\n")
  print(formula(theModel))
  
  if (length(intersect(class(theModel), c("glmerMod", "lmerMod", "lmerModLmerTest"))) > 0) {
    dropSummary <- data.frame(check.names = FALSE,
                              row.names = c("<none>", potentialDrops), 
                              Df = c("", sapply(res, function(x) attr(logLik(theModel), "df") - attr(logLik(x), "df"))),
                              #Deviance = round(c(deviance(theModel), sapply(X = res, FUN = deviance)), 4),
                              AIC = round(c(AIC(theModel), sapply(res, AIC)), 3),
                              LRT = c("", round(sapply(res, function(x) anova(x, theModel, test = "Chisq")$"Chisq"[2]), 5)),
                              `Pr(Chi)` = c("", round(sapply(res, function(x) anova(x, theModel, test = "Chisq")$"Pr(>Chisq)"[2]), 4)),
                              ` ` = c("", stars.pval(sapply(res, function(x) anova(x, theModel, test = "Chisq")$"Pr(>Chisq)"[2]))))
  } else {
    dropSummary <- data.frame(check.names = FALSE,
                              row.names = c("<none>", potentialDrops), 
                              df = c("", sapply(res, function(x) attr(logLik(theModel), "df") - attr(logLik(x), "df"))),
                              Deviance = round(c(deviance(theModel), sapply(X = res, FUN = deviance)), 4),
                              AIC = round(c(AIC(theModel), sapply(res, AIC)), 3),
                              LRT = c("", round(sapply(res, function(x) anova(x, theModel, test = "Chisq")$"Deviance"[2]), 5)),
                              `Pr(>Chi)` = c("", round(sapply(res, function(x) anova(x, theModel, test = "Chisq")$"Pr(>Chi)"[2]), 4)))   
  }
  
  if (length(intersect(class(theModel), c("glmerMod", "lmerMod", "lmerModLmerTest", "glm"))) > 0)
    print(dropSummary) else {
      
      print(cbind(dropSummary, Deviance = theDeviance)[, c("df", "Deviance", "AIC")])
    }
  
  # to-do: print method for lmer not correct (need Satterthwaite's via anova)
  cat("---\nSignif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1")
  
}

###
