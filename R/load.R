#' @export
Expr$methods(
  load = function() {
    filepaths <- sprintf("%s/%s/%s/%s.soft.gz",
                         organism, gpl,
                         valueType, gds)
    filepaths <- Filter(file.exists, filepaths)
    cat(sprintf("loading %s files\n", length(filepaths)))
    gdsData <<- plyr::llply(filepaths, softparser,
                      .parallel = T)
  }
)
