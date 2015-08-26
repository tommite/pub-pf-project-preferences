library(hitandrun)
library(plyr)

silta.constraints <- function() {
    n <- 6
    simplex <- simplexConstraints(n)

    w3w4.larger.than <- function(k) {
      cvec <- c(0, 0, -1, -1, 0, 0)
      cvec[k] <- 1
      list(constr=matrix(cvec, nrow=1), rhs=0, dir="<=")
    }
    w3w4.smaller.than <- function(k) {
      constr <- w3w4.larger.than(k)
      constr$constr <- -constr$constr
      constr
    }
    
    ## The criteria are: DPS Traffic Carry Width Salt Visual
    constr <- mergeConstraints(simplex,
                               ordinalConstraint(n, 1, 2),
                               ordinalConstraint(n, 2, 5),
                               ordinalConstraint(n, 2, 6),
                               w3w4.larger.than(5),
                               w3w4.larger.than(6),
                               w3w4.smaller.than(1))
    ## lower bounds of 0.02 on all criteria
    lb.constr <- do.call(mergeConstraints, llply(1:n, function(i) lowerBoundConstraint(n, i, 0.02)))
    constr <- mergeConstraints(constr, lb.constr)
    constr
}
