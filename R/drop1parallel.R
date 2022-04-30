library(doParallel)
library(MASS)
library(gtools)
library(lme4)

drop1parallel <- function (theModel, test = "Chisq", passData = NULL) { # only supports Chisq currently for all models (argument not used)
  
  if (!(class(theModel) %in% c("glmerMod", "clmm", "glmmTMB")))
    stop("You must supply a glmer, clmm, or glmmTMB model...")
  
  if (length(intersect(class(theModel), c("glmerMod", "lmerMod", "lmerModLmerTest"))) > 0) {
    
    if (is.null(passData))
      stop("You currently must pass data for a lme4 model...")
    
    if (class(theModel) == "glmerMod")
      originalData <- theModel@call[['data']]
    
    potentialDrops <- drop.scope(lme4:::getFixedFormula(formula(theModel)), formula("~ 1"))
    cat(paste("Dropping ", length(potentialDrops), " terms on ", getDoParWorkers(), " cores...\n\n", sep = ""))
    res <- foreach(theDrop = potentialDrops, .packages = attr(class(theModel), "package")) %dopar% update(theModel, as.formula(paste(".~.-", theDrop)), data = passData)
    
  } else if (class(theModel) %in% c("clmm", "glmmTMB")) {
    
    if (is.null(passData))
      stop("You currently must pass data for a clmm model...")
    else
      originalData <- passData
    
    if (class(theModel) == "clmm")
      requiredPackages <- "ordinal"
    else if (class(theModel) == "glmmTMB")
      requiredPackages <- "glmmTMB"
    
    potentialDrops <- drop.scope(lme4:::getFixedFormula(formula(theModel)), formula("~ 1"))
    cat(paste("Dropping ", length(potentialDrops), " terms on ", getDoParWorkers(), " cores...\n\n", sep = ""))
    res <- foreach(theDrop = potentialDrops, .packages = requiredPackages) %dopar% update(theModel, as.formula(paste(".~.-", theDrop)), data = originalData)
    
  } else {
    
    potentialDrops <- drop.scope(formula(theModel), formula("~ 1"))
    cat(paste("Dropping ", length(potentialDrops), " terms on ", getDoParWorkers(), " cores...\n\n", sep = ""))
    res <- foreach(theDrop = potentialDrops, .packages = attr(class(theModel), "package")) %dopar% update(theModel, as.formula(paste(".~.-", theDrop)))
    
  }
  
  if (length(intersect(class(theModel), c("glmerMod", "lmerMod", "lmerModLmerTest", "glm"))) > 0)
    for (each in 1:length(res)) {
      res[[each]]@call[['data']] <- originalData
    }
  
  if (length(intersect(class(theModel), c("glmmTMB"))) > 0)
    for (each in 1:length(res))
      res[[each]]$call[['data']] <- theModel$call$data
  #lapply(res, function(x) x@call[['data']])
  #anova(theModel, res[[1]])
  
  cat("Single term deletions\n\nModel:\n")
  print(formula(theModel))
  
  if (length(intersect(class(theModel), c("glmerMod", "lmerMod", "lmerModLmerTest", "glmmTMB"))) > 0) {
    dropSummary <- data.frame(check.names = FALSE,
                              row.names = c("<none>", potentialDrops), 
                              Df = c("", sapply(res, function(x) attr(logLik(theModel), "df") - attr(logLik(x), "df"))),
                              #Deviance = round(c(deviance(theModel), sapply(X = res, FUN = deviance)), 4),
                              AIC = round(c(AIC(theModel), sapply(res, AIC)), 3),
                              LRT = c("", round(sapply(res, function(x) anova(x, theModel, test = "Chisq")$"Chisq"[2]), 5)),
                              `Pr(Chi)` = c("", round(sapply(res, function(x) anova(x, theModel, test = "Chisq")$"Pr(>Chisq)"[2]), 4)),
                              ` ` = c("", gtools::stars.pval(sapply(res, function(x) anova(x, theModel, test = "Chisq")$"Pr(>Chisq)"[2])))) }
  else if (class(theModel) == "clmm") {
    
    dropSummary <- data.frame(check.names = FALSE,
                              row.names = c("<none>", potentialDrops), 
                              Df = c("", sapply(res, function(x) attr(logLik(theModel), "df") - attr(logLik(x), "df"))),
                              #Deviance = round(c(deviance(theModel), sapply(X = res, FUN = deviance)), 4),
                              AIC = round(c(AIC(theModel), sapply(res, AIC)), 3),
                              LRT = c("", round(sapply(res, function(x) ordinal:::anova.clm(x, theModel)["LR.stat"][2,1]), 3)),
                              `Pr(Chi)` = c("", round(sapply(res, function(x) ordinal:::anova.clm(x, theModel)["Pr(>Chisq)"][2,1]), 3)),
                              ` ` = c("", gtools::stars.pval(round(sapply(res, function(x) ordinal:::anova.clm(x, theModel)["Pr(>Chisq)"][2,1]), 4)))) }
  else {
    
    dropSummary <- data.frame(check.names = FALSE,
                              row.names = c("<none>", potentialDrops), 
                              df = c("", sapply(res, function(x) attr(logLik(theModel), "df") - attr(logLik(x), "df"))),
                              Deviance = round(c(deviance(theModel), sapply(X = res, FUN = deviance)), 4),
                              AIC = round(c(AIC(theModel), sapply(res, AIC)), 3),
                              LRT = c("", round(sapply(res, function(x) anova(x, theModel, test = "Chisq")$"Deviance"[2]), 5)),
                              `Pr(>Chi)` = c("", round(sapply(res, function(x) anova(x, theModel, test = "Chisq")$"Pr(>Chi)"[2]), 4)))   
  }
  
  if (length(intersect(class(theModel), c("glmerMod", "lmerMod", "lmerModLmerTest", "glm", "clmm", "glmmTMB"))) > 0)
    print(dropSummary) else {
      
      print(cbind(dropSummary, Deviance = theDeviance)[, c("df", "Deviance", "AIC")])
    }
  
  # to-do: print method for lmer not correct (need Satterthwaite's via anova)
  cat("---\nSignif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n\n")
  
  return(dropSummary)
  
}

###

stepAICparallel <- function (startMod, theThreshold = 0, passData = NULL, ignoreTerms = NULL, evaluateTerms = NULL) {
  
  start.time <- Sys.time()
  
  startPoint <- drop1parallel(startMod, test = "Chisq", passData)
  
  startFrame <- data.frame(term = rownames(startPoint), pVal = startPoint[, grep("Pr", names(startPoint))], AIC = startPoint[, grep("AIC", names(startPoint))])[-1,]
  rownames(startFrame) <- NULL
  
  potentialDrops <- drop.scope(lme4:::getFixedFormula(formula(startMod)), formula("~ 1"))
  
  if (!is.null(ignoreTerms)) {
    if (sum(ignoreTerms %in% potentialDrops) != length(ignoreTerms))
      stop("One or more ignore terms not found in the model...")
    cat(paste("\nIgnoring specified terms:", paste(ignoreTerms, collapse = ", "), "\n"))
    
    startFrame <- subset(startFrame, !(term %in% ignoreTerms))
  }
  
  if (!is.null(evaluateTerms)) {
    if (sum(evaluateTerms %in% potentialDrops) != length(evaluateTerms))
      stop("One or more evaluate terms not found in the model...")
    cat(paste("\nEvaluating specified terms:", paste(evaluateTerms, collapse = ", "), "\n"))
    
    startFrame <- subset(startFrame, term %in% evaluateTerms)
  }
  
  if (length(intersect(ignoreTerms, evaluateTerms)) > 0)
    stop("Ignored and evaluated terms should not intersect...")
  
  # initialize
  stepDownMod <- startMod
  finalStepLRT <- startPoint
  dropSummary <- data.frame(dropTerm = NA, df = NA, AIC = NA)
  
  loopStore <- list(startMod)
  dropStore <- list(NA)
  
  # to-do: need to warn/error if REML = FALSE not set for lmer
  # to-do: this won't work for lmer until drop1parallel returns p values
  
  theThreshold <- 0
  
  if (AIC(startMod) - min(startFrame$AIC) > theThreshold) {
    
    paste("Next Drop:", subset(startFrame, AIC == min(AIC))$term)
    
    repeat {
      
      if ( AIC(stepDownMod) - min(startFrame$AIC) < theThreshold ) {
        cat("\nNothing left to drop! Stopping...\n\n")
        break
      }
      
      fieldToRemove <- as.character(subset(startFrame, AIC == min(AIC))$term)
      
      cat(paste("Now Dropping: ", fieldToRemove, "\n", sep=""))
      
      stepDownMod <- update(stepDownMod, as.formula(paste(".~.-", fieldToRemove)) )
      
      loopStore[[length(loopStore) + 1]] <- stepDownMod
      dropStore[[length(dropStore) + 1]] <- fieldToRemove
      
      startPoint <- drop1parallel(stepDownMod, test = "Chisq", passData)
      
      startFrame <- data.frame(term = rownames(startPoint), pVal = startPoint[, grep("Pr", names(startPoint))], AIC = startPoint[, grep("AIC", names(startPoint))])[-1,]
      rownames(startFrame) <- NULL
      
      if (!is.null(ignoreTerms))
        startFrame <- subset(startFrame, !(term %in% ignoreTerms))
      
      if (!is.null(evaluateTerms))
        startFrame <- subset(startFrame, term %in% evaluateTerms)
      
      startFrame[order(startFrame$pVal),]
      
    }
    
    AIC(startMod, stepDownMod)
    
    theDFs <- list()
    for (i in 1:length(loopStore))
      theDFs[i] <- attr(logLik(loopStore[[i]]), "df")
    
    (finalStepLRT <- startPoint)
    
    (dropSummary <- data.frame(dropTerm = unlist(dropStore), df = unlist(theDFs), AIC = sapply(X = loopStore, FUN = AIC)))
    
  } else cat("No Potential Drops! Aborting...")
  
  return(list(stepDownMod, finalStepLRT, dropSummary)) # should name list entries
}
