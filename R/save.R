#' @export
Expr$methods(
  save = function(file) {
    'Save the current object on the file
    in R external object format.
    '
    cat(sprintf("saving object\n"))
    base::saveRDS(object = .self,
                  file = paste(file,"rds"
                               ,sep = "."))}
)
