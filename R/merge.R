#' @export
Expr$methods(
  merge = function(fbind = FALSE) {
    if(!fbind){
      cat(sprintf("merging %s tables\n", length(gdsData)))
      gdsData <<- Reduce(function(x,y)
        x[y], gdsData)
    }else if(fbind == TRUE) {
      if(class(gdsData) == "list" )
        cat(sprintf("merging %s tables\n", length(gdsData)))
        gdsData <<- do.call(cbind,gdsData)
        gdsData <<- gdsData[,unique(names(gdsData)),with=F]
      }
  }
)
