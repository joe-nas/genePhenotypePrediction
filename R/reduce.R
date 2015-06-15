#'@export
Expr$methods(
  reduce = function(){
    gdssvd <- svd(gdsData[,3:ncol(gdsData), with = F], nv = ncol(gdsData)-2, nu = 0)
    var_explained <<- data.frame(PCs = 1:length(gdssvd$d),
                                 Var = cumsum(gdssvd$d/sum(gdssvd$d)))
    slope <- diff(var_explained[,"Var"])/diff(var_explained[,"PCs"])
    n_pcs <<- max(which(slope > quantile(slope,.85)))
    red <- quote((as.matrix(gdsData[,3:ncol(gdsData),with=F]) %*% t(gdssvd$v))[,1:n_pcs])
    gdsRed <<- data.table(gdsData[,1:2, with=F],eval(red))
  }
)
