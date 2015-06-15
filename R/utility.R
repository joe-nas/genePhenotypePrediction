#' @export
softparser <- function(filepath){
  require(data.table)
  tmp <- readLines(filepath)
  tmp <- tmp[(grep("^!dataset_table_begin", tmp)+1):(grep("^!dataset_table_end", tmp)-1)]
  tmpcolumns <- unlist(strsplit(tmp[1],"\t"))
  gsmcolumns <- sum(grepl("GSM",tmpcolumns))
  dt <- fread(paste0(gsub("null","NA",tmp), collapse = "\n"),
              colClasses = c(rep("character",length(tmpcolumns)-gsmcolumns),
                             rep("numeric",gsmcolumns)),
              na.strings = 'null', sep="\t", header = T,
              integer64 = "double",
              verbose = F)[, .SD,, by=c("ID_REF","IDENTIFIER")]
  dt <- setkey(dt, IDENTIFIER, ID_REF)[]
  return(dt[order(IDENTIFIER)][])}
