#'@export
Expr$methods(
  reduce = function(slope = .80){
    if(is.data.table(gdssvd) == FALSE){
      cat(sprintf("svd\n"))
      svd.res <- svd(gdsData[,3:ncol(gdsData), with = F],
                     nv = ncol(gdsData)-2, nu = 0)
      rv <<- svd.res$v
      sv <<- svd.res$d
      cat(sprintf("projecting\n"))
      gdssvd <<- data.table(as.matrix(gdsData[,3:ncol(gdsData),with=F]) %*% t(rv))
      }

    var_explained <<- data.frame(PCs = 1:length(sv),
                                 Var = cumsum(sv/sum(sv)))
    sl <- diff(var_explained[,"Var"])/diff(var_explained[,"PCs"])

    n_pcs <<- max(which(sl > quantile(sl,slope)))

    gdsRed <<- data.table(gdsData[,1:2, with=F],gdssvd[,1:n_pcs,with=F])
    invisible(gdsRed)
  }
)
