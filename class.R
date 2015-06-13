library(GEOmetadb)
library(GEOquery)
library(plyr)
library(data.table)
library(caret)
library(doMC)
registerDoMC(10)

con <- dbConnect(SQLite(),'GEOmetadb.sqlite')
gdsquery <- dbGetQuery(con,'SELECT sample_organism, gds, sample_count, gpl, value_type from gds WHERE (sample_organism="Homo sapiens" OR sample_organism="Mus musculus") AND (value_type="count" OR value_type="transformed count")')
dbDisconnect(con)

# table(dbGetQuery(con,'SELECT * from gds WHERE (sample_organism="Homo sapiens" OR sample_organism="Mus musculus")')$value_type)


gdsquery <- data.frame(lapply(gdsquery,function(x) gsub(" ","_",x)))
gdsquery$sample_count <- as.integer(gdsquery$sample_count)
gds <- Filter(function(x){sum(x$sample_count) >= 50},split(gdsquery, list(gdsquery$sample_organism, gdsquery$gpl, gdsquery$value_type),))

soft2dt <- function(filepath){
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
              verbose = F)[, .SD,, c("ID_REF","IDENTIFIER")]
  return(dt)
}


setClassUnion("Expr.gds.dataTypes", c("list", "data.frame"))
Expr <- setRefClass(Class = "Expr",
                    fields = list(
                      organism = "character",
                      valueType = "character", 
                      gpl = "character", 
                      gds = "character",
                      sampleCounts = "numeric", 
                      gdsData = "Expr.gds.dataTypes",
                      gdsRed = "data.table"),
                    
                    
                    methods = list(
                      download = function(){
                        return()
                        },
                      load = function(){
                        filepaths <- sprintf("%s/%s/%s/%s.soft.gz", 
                                             organism, gpl, 
                                             valueType, gds)
                        filepaths <- Filter(file.exists, filepaths)
                        gdsData <<- llply(filepaths, soft2dt,
                                          .parallel = T)
                        },
                      merge = function(){
                        gdsData <<- Reduce(function(x,y) x[y], gdsData)
                      },
                      impute = function(){
                        if(is.data.table(gdsData)){
                          preproc <- preProcess(gdsData[,3:ncol(gdsData),with=F], 
                                                method = c("center","scale","medianImpute"))
                          gdsRed <<- gdsData[,3:ncol(gdsData) := data.table(predict(preproc,gdsData[,3:ncol(gdsData),with=F])), with=F]
                        }else if(is.list(gdsData)){
                          cat("use $merge() before impute()")
                        }},
                      reduce = function(n_pcs = FALSE){
                        gdssvd <- svd(gdsData[,3:ncol(gdsData), with = F], nv = ncol(gdsData)-2, nu = 0)
                        var_explained <- data.table(PCs = 1:length(gdssvd$d), Var = cumsum(gdssvd$d/sum(gdssvd$d)))
                        if(!n_pcs){
                          slope <- diff(var_explained[,Var])/diff(var_explained[,PCs])
                          n_pcs <- max(which(slope > quantile(slope,.90)))
                        }
                        gdsRed <<- data.table(gdsData[,1:2, with=F],data.table((as.matrix(gdsData[,3:ncol(gdsData),with=F]) %*% t(gdssvd$v))[,1:n_pcs]))
                      }
                      )
)

initializeExpr <- function(gds = "list"){
  lapply(gds, function(o){ 
    Expr( organism = gsub(" ","_",as.character(unique(o$sample_organism))),
          valueType = gsub(" ","_",as.character(unique(o$value_type))),
          gpl = as.character(unique(o$gpl)),
          gds = as.character(unique(o$gds)),
          sampleCounts = as.numeric(o$sample_counts))
  })
}

test <- initializeExpr(gds)
l_ply(1:10,function(i){
  test[[i]]$load()
  test[[i]]$merge()
  test[[i]]$impute()
  test[[i]]$reduce()
})
