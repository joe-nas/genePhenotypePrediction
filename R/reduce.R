#'@export
Expr$methods(
  reduce = function(slope = .80){
    exprSvd <<- new("ExprSVD",
                    "organism" = organism,
                    "gpl" = gpl,
                    "valueType" = valueType)

    if(all(complete.cases(gdsData))){
      dat <- quote(gdsData[,3:ncol(gdsData),with=F])
    }else{
      imp <- .self$impute()
      dat <- quote(rbindlist(list(gdsData, imp),
                       use.names = T)[,3:ncol(gdsData),with=F])
    }
    dosvd <- quote(
      svd(eval(dat),
      nv = ncol(gdsData)-2, nu = 0))


    with(eval(dosvd), {
      exprSvd$d <<- d
      cat(sprintf("svd\n"))
      dat <- as.matrix(eval(dat))
      cat(sprintf("projecting\n"))
      exprSvd$XV <<- Matrix::tcrossprod(
        # Matrix::Matrix(as.matrix(gdsData[,3:ncol(gdsData),with=F])),
        Matrix::Matrix(dat),
        Matrix::Matrix(v))
      })


    exprSvd$Var <<- data.frame(
      PCs = 1:length(exprSvd$d),
      Var = cumsum(exprSvd$d/sum(exprSvd$d)))


    sl <- diff(
      exprSvd$Var[,"Var"])/diff(exprSvd$Var[,"PCs"])

    n_pcs <<- max(
      which(sl > quantile(sl,slope)))

    invisible(exprSvd)
  }
)
