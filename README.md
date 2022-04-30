# drop1parallel
A parallel version of the drop1 (and stepAIC) function for R

### Installation
```R
library(devtools)
install_github("https://github.com/davedgd/drop1parallel")
```

### Usage
```R
library(lme4)
library(doParallel)
library(HSAUR3)
library(drop1parallel)

cl <- makePSOCKcluster(detectCores(logical = FALSE))
registerDoParallel(cl)

# basic example from ?lme4::glmer
gm2 <- glmer(outcome~ treatment + poly(visit, 2) +(1|patientID),
               data=toenail,
               family=binomial,nAGQ=20)

stepAICparallel(gm2, passData = toenail)
```

### To-Do
- code cleanup
- add an option for automated polynomial reduction (e.g., cubic -> quadratic -> linear)
- report convergence or model fit issues during stepwise reduction
