#' @importClassesFrom  data.table data.table
#' @importClassesFrom  Matrix dgeMatrix
#' @export Expr.gds.dataTypes Expr.gdssvd.dataTypes
#' @exportClass Expr



Expr.gds.dataTypes <- setClassUnion("Expr.gds.dataTypes",
                                    c("list", "data.frame"))
Expr.gdssvd.dataTypes <- setClassUnion("Expr.gdssvd.dataTypes",
                                       c("logical", "data.table",
                                         "dgeMatrix"))

Expr <- setRefClass(
  Class = "Expr",
  fields = list(
    organism = "character",
    valueType = "character",
    gpl = "character",
    gds = "character",
    sampleCounts = "integer",
    gdsData = "Expr.gds.dataTypes",
    gdsMissing = "data.table",
    thresh = "ANY",
    exprSvd = "ExprSVD",
    n_pcs = "integer",
    fit  = "ANY",
    prediction = "ANY"
  )#,
#   methods = list(
#     initialize = function(...){
#       gdssvd <<- FALSE
#       callSuper(...)
#       }
#     )
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
