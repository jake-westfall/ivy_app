
semPlotModel_MxRAMModel <- function (object) 
{
  varNames <- object@manifestVars
  factNames <- object@latentVars
  Dirpaths <- which(t(object@matrices$A@free | object@matrices$A@values != 
                        0), arr.ind = TRUE)
  DirpathsFixed <- !t(object@matrices$A@free)[Dirpaths]
  DirpathsValues <- t(object@matrices$A@values)[Dirpaths]
  DirpathsLabels <- t(object@matrices$A@labels)[Dirpaths]
  Sympaths <- which(t(object@matrices$S@free | object@matrices$S@values != 
                        0) & upper.tri(object@matrices$S@values, diag = TRUE), 
                    arr.ind = TRUE)
  SympathsFixed <- !t(object@matrices$S@free)[Sympaths]
  SympathsValues <- t(object@matrices$S@values)[Sympaths]
  SympathsLabels <- t(object@matrices$A@labels)[Sympaths]
  if (!is.null(object@matrices$M)) {
    Means <- which(object@matrices$M@free | object@matrices$M@values != 
                     0)
    MeansFixed <- !object@matrices$M@free[Means]
    MeansValues <- object@matrices$M@values[Means]
    MeansLabels <- object@matrices$M@labels[Means]
  }
  else {
    Means <- numeric(0)
    MeansFixed <- logical(0)
    MeansValues <- numeric(0)
    MeansLabels <- character(0)
  }
  if (!length(object@output) == 0) {
    standObj <- standardizeRam(object, "model")
    DirpathsValuesStd <- t(standObj@matrices$A@values)[Dirpaths]
    SympathsValuesStd <- t(standObj@matrices$S@values)[Sympaths]
    if (!is.null(standObj@matrices$M)) {
      MeansValuesStd <- standObj@matrices$S@values[Means]
    }
    else {
      MeansValuesStd <- numeric(0)
    }
  }
  else {
    DirpathsValuesStd <- rep(NA, nrow(Dirpaths))
    SympathsValuesStd <- rep(NA, nrow(Sympaths))
    MeansValuesStd <- rep(NA, length(Means))
  }
  Vars <- data.frame(name = c(varNames, factNames), manifest = c(varNames, 
                                                                 factNames) %in% varNames, exogenous = NA, stringsAsFactors = FALSE)
  Pars <- data.frame(label = c(DirpathsLabels, SympathsLabels, 
                               MeansLabels), lhs = c(Vars$name[c(Dirpaths[, 1], Sympaths[, 
                                                                                         1])], rep("", length(Means))), edge = c(rep("->", nrow(Dirpaths)), 
                                                                                                                                 rep("<->", nrow(Sympaths)), rep("int", length(Means))), 
                     rhs = Vars$name[c(Dirpaths[, 2], Sympaths[, 2], Means)], 
                     est = c(DirpathsValues, SympathsValues, MeansValues), 
                     std = c(DirpathsValuesStd, SympathsValuesStd, MeansValuesStd), 
                     group = object@name, fixed = c(DirpathsFixed, SympathsFixed, 
                                                    MeansFixed), par = 0, stringsAsFactors = FALSE)
  Pars$par[is.na(Pars$label)] <- seq_len(sum(is.na(Pars$label)))
  for (lbl in unique(Pars$label[!is.na(Pars$label)])) {
    Pars$par[Pars$label == lbl] <- max(Pars$par) + 1
  }
  Pars$label[is.na(Pars$label)] <- ""
  semModel <- new("semPlotModel")
  semModel@Pars <- Pars
  semModel@Vars <- Vars
  semModel@Computed <- !length(object@output) == 0
  semModel@Original <- list(object)
  if (!is.null(object@data)) {
    if (object@data@type == "cov") {
      semModel@ObsCovs <- list(object@data@observed)
    }
    else if (object@data@type == "raw") {
      semModel@ObsCovs <- list(cov(object@data@observed))
    }
    else {
      semModel@ObsCovs <- list(NULL)
    }
  }
  else {
    semModel@ObsCovs <- list(NULL)
  }
  semModel@ImpCovs <- list(object$fitfunction$info$expCov)
  return(semModel)
}

standardizeRam <- function (model, return = "parameters", Amat = NA, Smat = NA, 
          Mmat = NA) 
{
  if (!(return == "parameters" | return == "matrices" | return == 
          "model")) 
    stop("Invalid 'return' parameter. What do you want from me?")
  obj <- class(model$expectation)[1]
  suppliedNames <- !is.na(Amat) & !is.na(Smat)
  cA <- is.character(Amat)
  cS <- is.character(Smat)
  cM <- is.character(Mmat)
  if (obj != "MxExpectationRAM" & (!cA)) 
    stop("I need either MxExpectationRAM or the names of the A and S matrices.")
  output <- model@output
  if (is.null(output)) 
    stop("Provided model has no objective function, and thus no output. I can only standardize models that have been run!")
  if (length(output) < 1) 
    stop("Provided model has no output. I can only standardize models that have been run!")
  if (cA) {
    nA <- Amat
  }
  else {
    nA <- model$expectation$A
  }
  if (cS) {
    nS <- Smat
  }
  else {
    nS <- model$expectation$S
  }
  if (cM) {
    nM <- Mmat
  }
  else {
    nM <- model$expectation$M
  }
  A <- model[[nA]]
  S <- model[[nS]]
  d <- dim(S@values)[1]
  I <- diag(d)
  IA <- solve(I - A@values)
  expCov <- IA %*% S@values %*% t(IA)
  invSDs <- 1/sqrt(diag(expCov))
  names(invSDs) <- as.character(1:length(invSDs))
  if (!is.null(dimnames(A@values))) {
    names(invSDs) <- as.vector(dimnames(S@values)[[2]])
  }
  diag(I) <- invSDs
  stdA <- I %*% A@values %*% solve(I)
  stdS <- I %*% S@values %*% I
  model[[nA]]@values[, ] <- stdA
  model[[nS]]@values[, ] <- stdS
  if (!is.na(nM)) {
    model[[nM]]@values[, ] <- rep(0, length(invSDs))
  }
  if (return == "model") 
    return(model)
  matrices <- list(model[[nA]], model[[nS]])
  names(matrices) <- c("A", "S")
  if (return == "matrices") 
    return(matrices)
  p <- summary(model)$parameters
  p <- p[(p[, 2] == nA) | (p[, 2] == nS), ]
  rescale <- invSDs[p$row] * 1/invSDs[p$col]
  rescaleS <- invSDs[p$row] * invSDs[p$col]
  rescale[p$matrix == "S"] <- rescaleS[p$matrix == "S"]
  p[, 5] <- p[, 5] * rescale
  p[, 6] <- p[, 6] * rescale
  names(p)[5:6] <- c("Std. Estimate", "Std.Std.Error")
  return(p)
}