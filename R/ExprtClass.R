#' @importClassesFrom  data.table data.table
#' @export Expr.gds.dataTypes
#' @exportClass Expr



Expr.gds.dataTypes <-  setClassUnion("Expr.gds.dataTypes",
                                     c("list", "data.frame"))

Expr <- setRefClass(
  Class = "Expr",
  fields = list(
    organism = "character",
    valueType = "character",
    gpl = "character",
    gds = "character",
    sampleCounts = "integer",
    gdsData = "Expr.gds.dataTypes",
    gdsRed = "data.table",
    n_pcs = "integer",
    var_explained = "data.frame",
    #prediction = "data.frame"
    prediction = "ANY"
  )
)
#
# ResultClass <- setRefClass(
#   Class = "ResultClass",
#   fields = list(
#     var_explained = "data.frame",
#     n_pcs = "interger",
#     fit = "list",
#     prediction = "data.frame",
#   )
# )
