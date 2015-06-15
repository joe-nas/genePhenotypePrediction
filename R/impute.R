#'@export
Expr$methods(
  impute = function(method,maxmissing = 0.5){
    if(is.data.table(gdsData)){
      gdsRed <<- gdsData[!rowSums(is.na(gdsData)) >ncol(gdsData)*maxmissing,]
      gdsData <<- gdsData[!rowSums(is.na(gdsData)) >ncol(gdsData)*maxmissing,]
#      missing_ij <- which(is.na(test[[1]]$gdsData[[1]]),arr.ind = T)
#       for(i in missing_ij){
#         set(gdsData,i=mising_ij[i,1],j=mising_ij[i,2],value = j=mising_ij[i,3])
#       }
      preproc <- caret::preProcess(gdsRed[,3:ncol(gdsRed),with=F],
                                   method = c("center","scale","medianImpute"))
      cat(sprintf("imputing\n"))
      imputed <- quote(data.table(predict(preproc, gdsRed[,3:ncol(gdsRed),with=F])))
      gdsRed <<- gdsData[,3:ncol(gdsRed) := eval(imputed) , with=F]
    }else if(is.list(gdsData)){
      cat("use $merge() before impute()")
    }}
)
