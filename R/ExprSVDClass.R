#' @importClassesFrom  Matrix dgeMatrix

ExprSVD <- setRefClass(
  Class = "ExprSVD",
  fields = list(
    organism = "character",
    gpl = "character",
    valueType = "character",
    XV = "dgeMatrix",
    d = "numeric",
    Var = "data.frame"
    ),
  methods = list(
    save = function(){
      cat(sprintf("saving: %s %s %s\n",
              organism,gpl,valueType))
      base::saveRDS(
        object = .self,
        file = sprintf("%s_%s_%s.rds",
                       organism,gpl,valueType)
        )
    })
  )


ExprPrediction <- setRefClass(
  Class = "ExprPrediction",
  fields = list(
    splitPrediction = "list",
    combinedPrediction = "data.table"
  )
)
