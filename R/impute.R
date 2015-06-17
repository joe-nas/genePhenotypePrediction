#'@export
Expr$methods(
  impute = function(method, maxmissing = 0.5){
    if(is.data.table(gdsData)){
      if(is.null(gdsMissing)){
        complete <- complete.cases(gdsData)
        if(all(complete)){
          return(invisible(gdsData))
        }
        gdsMissing <<- gdsData[which(!complete),]
        gdsData <<- gdsData[which(complete),]
        }

      dtdim <- dim(gdsMissing)
      thresh <<- gdsMissing[rowSums(is.na(gdsMissing)) <= (dtdim[2]-2)*maxmissing,]

      if(nrow(thresh) > 0){
        miss <- is.na(thresh[,3:dtdim[2],with=F])
        rowmed <- suppressWarnings(
          thresh[,median(as.numeric(.SD),na.rm=T),by=1:nrow(thresh)][,V1])
        values <- data.table(
          as.matrix(gdata::NAToUnknown(thresh[,3:dtdim[2],with=F],0)) + rowmed * miss)
#         imputed <- quote(predict(
#           caret::preProcess(thresh[,3:dtdim[2], with=F],
#                             method = c("center", "scale", "medianImpute")),
#           newdata = thresh[,3:dtdim[2], with=F]))
        cat(sprintf("imputing\n"))
        return(invisible(thresh[,(3:dtdim[2]):=values][,.SD]))
      }else{
        return(invisible(gdsData))
      }
    }else if(is.list(gdsData)){
      cat("use $merge() before impute()")
    }}
)
