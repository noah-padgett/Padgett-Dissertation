# categorical omega


.catOmega <- function(dat) {

  if(!requireNamespace("lavaan", quietly = TRUE)) stop("The package 'lavaan' is needed; please install the package and try again.")
  if(!requireNamespace("mnormt", quietly = TRUE)) stop("The package 'mnormt' is needed; please install the package and try again.")

  q <- ncol(dat)
  for(i in 1:q) dat[,i] <- ordered(dat[,i])
  varnames <- paste0("y", 1:q)
  colnames(dat) <- varnames
  loadingName <- paste("a", 1:q, sep = "")
  errorName <- paste("b", 1:q, sep = "")
  model <- "f1 =~ NA*y1 + "
  loadingLine <- paste(paste(loadingName, "*", varnames, sep = ""), collapse = " + ")
  factorLine <- "f1 ~~ 1*f1\n"
  model <- paste(model, loadingLine, "\n", factorLine)
  error <- try(fit <- lavaan::cfa(model, data = dat, se = "none", ordered = varnames), silent = TRUE)
  converged <- FALSE
  if(!is(error, "try-error") && fit@Fit@converged) converged <- TRUE
  reliab <- NA
  if(converged) {
    param <- lavaan::inspect(fit, "coef")
    ly <- param[["lambda"]]
    ps <- param[["psi"]]
    truevar <- ly%*%ps%*%t(ly)
    threshold <- .getThreshold(fit)[[1]]
    denom <- .polycorLavaan(fit, dat)[varnames, varnames]
    invstdvar <- 1 / sqrt(diag(fit@Fit@Sigma.hat[[1]]))
    polyr <- diag(invstdvar) %*% truevar %*% diag(invstdvar)
    print(polyr)
    sumnum <- 0
    addden <- 0
    for(j in 1:q) {
      for(jp in 1:q) {
        sumprobn2 <- 0
        addprobn2 <- 0
        t1 <- threshold[[j]]
        t2 <- threshold[[jp]]
        for(c in 1:length(t1)) {
          for(cp in 1:length(t2)) {
            sumprobn2 <- sumprobn2 + .p2(t1[c], t2[cp], polyr[j, jp])
            addprobn2 <- addprobn2 + .p2(t1[c], t2[cp], denom[j, jp])
          }
        }
        sumprobn1 <- sum(pnorm(t1))
        sumprobn1p <- sum(pnorm(t2))
        sumnum <- sumnum + (sumprobn2 - sumprobn1 * sumprobn1p)
        addden <- addden + (addprobn2 - sumprobn1 * sumprobn1p)
      }
    }
    reliab <- sumnum / addden
  }
  reliab
}

.p2 <- function(t1, t2, r) {
  mnormt::pmnorm(c(t1, t2), c(0,0), matrix(c(1, r, r, 1), 2, 2))
}


.polycorLavaan <- function(object, data) {
  ngroups <- object@Data@ngroups
  coef <- lavaan::inspect(object, "coef")
  targettaunames <- NULL
  if(ngroups == 1) {
    targettaunames <- rownames(coef$tau)
  } else {
    targettaunames <- rownames(coef[[1]]$tau)
  }
  barpos <- sapply(strsplit(targettaunames, ""), function(x) which(x == "|"))
  varnames <- unique(apply(data.frame(targettaunames, barpos - 1), 1, function(x) substr(x[1], 1, x[2])))
  script <- ""
  for(i in 2:length(varnames)) {
    temp <- paste0(varnames[1:(i - 1)], collapse = " + ")
    temp <- paste0(varnames[i], "~~", temp, "\n")
    script <- paste(script, temp)
  }
  suppressWarnings(newobject <- .refit(script, data, varnames, object))
  if(ngroups == 1) {
    return(lavaan::inspect(newobject, "coef")$theta)
  } else {
    return(lapply(lavaan::inspect(newobject, "coef"), "[[", "theta"))
  }
}

.getThreshold <- function(object) {
  ngroups <- object@Data@ngroups
  coef <- lavaan::inspect(object, "coef")
  result <- NULL
  if(ngroups == 1) {
    targettaunames <- rownames(coef$tau)
    # This simply identifies where the "|" for separating the variable names from the thresholds.
    barpos <- sapply(strsplit(targettaunames, ""), function(x) which(x == "|"))
    varthres <- apply(data.frame(targettaunames, barpos-1), 1, function(x) substr(x[1], 1, x[2]))
    result <- list(split(coef$tau, factor(varthres, levels=unique(varthres))))
  } else {
    result <- list()
    for(g in 1:ngroups) {
      stop(print("Developmental; please contact maintainer!"))
      #targettaunames <- rownames(coef[[g]]$tau)
      #barpos <- sapply(strsplit(targettaunames, ""), function(x) which(x == "|"))
      #varthres <- apply(data.frame(targettaunames, barpos-1), 1, function(x) substr(x[1], 1, x[2]))
      #print(varthres)
      #result[[g]] <- split(coef[[g]]$tau, factor(varthres, levels=unique(varthres)))
    }
  }
  return(result)
}

.refit <- function(pt, data, vnames, object) {
  previousCall <- object@call
  args <- as.list(previousCall[-1])
  args$model <- pt
  args$data <- data
  args$ordered <- vnames
  funcall <- as.character(previousCall[[1]])
  tempfit <- do.call(funcall[length(funcall)], args)
}

.catOmega(mydata[,1:10])


library(MBESS)

fit.con <- MBESS::ci.reliability(mydata[,1:10], type="omega")
fit.cat <- MBESS::ci.reliability(mydata[,1:10], type="cat",interval.type = "bca",B = 100)
fit.con$est
