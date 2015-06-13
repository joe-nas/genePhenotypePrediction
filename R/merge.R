#' @export
Expr$methods(
  merge = function() {
    gdsData <<- Reduce(function(x,y)
      x[y], gdsData)
  }
)
