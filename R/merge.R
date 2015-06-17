#' @export
Expr$methods(
  merge = function(fbind = FALSE) {
    if(class(gdsData) == "list" ){
      if(!fbind){
        cat(sprintf("merging %s tables\n", length(gdsData)))
        gdsData <<- Reduce(function(x,y)
          x[y], gdsData)
      }else if(fbind == TRUE) {
        cat(sprintf("merging %s tables\n", length(gdsData)))
        gdsData <<- do.call(base::cbind, gdsData)
        gdsData <<- gdsData[,unique(names(gdsData)),with=F]
      }
    }else{
      sprintf("nothing to merge")
    }
  }
)
