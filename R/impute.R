#'@export
Expr$methods(
  impute = function(){
    if(is.data.table(gdsData)){
      preproc <- caret::preProcess(gdsData[,3:ncol(gdsData),with=F],
                                   method = c("center","scale","medianImpute"))
      imputed <- quote(data.table(predict(preproc, gdsData[,3:ncol(gdsData),with=F])))
      gdsRed <<- gdsData[,3:ncol(gdsData) := eval(imputed) , with=F]
    }else if(is.list(gdsData)){
      cat("use $merge() before impute()")
    }}
)
