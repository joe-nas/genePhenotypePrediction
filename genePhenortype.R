library(GEOmetadb)
library(GEOquery)
library(plyr)
# download database
# if(!file.exists('GEOmetadb.sqlite')) getSQLiteFile()



# GSM sample, GPL platform, GSE series, GDS curated
# dbListFields(con,'gpl')


# gsm, gpl, gds
# menschundmaus <- dbGetQuery(con,'SELECT DISTINCT gpl.bioc_package, gpl.organism, gpl.gpl ,gse.gse FROM gse_gpl INNER JOIN gpl ON gse_gpl.gpl=gpl.gpl INNER JOIN gse ON gse_gpl.gse=gse.gse  WHERE bioc_package IS NOT NULL AND (organism="Homo sapiens" OR organism="Mus musculus")')
# menschundmaus_l <- Filter(nrow,split(menschundmaus, list(menschundmaus$organism,menschundmaus$gpl)))
con <- dbConnect(SQLite(),'GEOmetadb.sqlite')
gdsquery <- dbGetQuery(con,'SELECT sample_organism, gds, sample_count, gpl, value_type from gds WHERE (sample_organism="Homo sapiens" OR sample_organism="Mus musculus") AND (value_type="count" OR value_type="transformed count")')
dbDisconnect(con)

# table(dbGetQuery(con,'SELECT * from gds WHERE (sample_organism="Homo sapiens" OR sample_organism="Mus musculus")')$value_type)


gdsquery <- data.frame(lapply(gdsquery,function(x) gsub(" ","_",x)))
gdsquery$sample_count <- as.integer(gdsquery$sample_count)
gds <- Filter(function(x){sum(x$sample_count) >= 50},split(gdsquery, list(gdsquery$sample_organism, gdsquery$gpl, gdsquery$value_type),))



downLoadGDS <- function(gds, download_dest = ".", method = "dl" ){
  gds_l <- lapply(gds, function(z) gsub(" ","_", z))
  path <- sprintf("%s/%s/%s/%s", download_dest, gds_l$sample_organism, gds_l$gpl, gds_l$value_type)
  cat(path,"\n")
  dir.create(path, recursive = TRUE)
  out = NA
  if(method == "d"){
    try({out<-getGEOfile(gds$gds,destdir = path)}, silent = F)
  }
  else if(method == "l"){
    try({out<-getGEO(gds$gds,GSEMatrix = T,filename = sprintf("%s/%s.soft.gz", path, gds$gds)) }, silent = F)
  }  
  else if(method == "dl"){
    try({out<-getGEO(GEO = gds$gds, GSEMatrix = T, destdir = path, filename = sprintf("%s.soft.gz", gds$gds)) }, silent = F)
  }
  out
}



# downloaden
library(doMC)
registerDoMC(4)
expression_data <- llply(seq_along(gds), function(i){
  cat(gds[[i]]$gpl,":")
  alply(gds[[i]], 1, function(y){
    cat(y$sample_count," ")
    downLoadGDS(y,".","d")
}, .inform = T,.parallel = T)})

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
              verbose = T)[, .SD,, c("ID_REF","IDENTIFIER")]
  return(dt)
}

softpath <- function(workingdir,gdsline){
  sprintf("%s/%s/%s/%s/%s.soft.gz", workingdir, gdsline$sample_organism, gdsline$gpl, gdsline$value_type, gdsline$gds)
}

gpldtpath <- function(workingdir, gdsline){
  sprintf("%s/%s/%s_%s_%s_table.csv", workingdir, gdsline$sample_organism, 
          gdsline$sample_organism, gdsline$gpl, gdsline$value_type)
}

savedt <- function(gpldtpath, dt){
  write.table(dt ,file = gpldtpath, col.names = T, row.names = F, sep = ",")
}

gdspaths <- Map(function(x)softpath(".",x),gds)
fexists <- Map(file.exists,gdspaths)
gdspaths <- Filter(length,Map(function(x) gdspaths[[x]][fexists[[x]]],seq_along(gdspaths)))

registerDoMC(3)
gdspaths <- gdspaths[107:111]
res <- Map(function(i) Reduce(function(x,y) x[y],Map(soft2dt,gdspaths[[i]])),1:length(gdspaths))
# res.copy <- res
# laply(res,dim)
# res <- res[[2]]

#impute.knn(as.matrix(res[,(.SD),.SDcols = grep("GSM",names(res))]))
library(caret)
library(lattice)

preprocessData <- function(dt){
  dt_dim <- dim(dt)
  PreProcValues <- preProcess(dt[,(.SD),.SDcols = grep("GSM",names(dt))], method = c("center","scale","medianImpute"))
  transformedData <- data.table(predict(PreProcValues, dt[,(.SD),.SDcols = grep("GSM",names(dt))]))
  dat_svd <- svd(transformedData, nu = nrow(transformedData),nv = ncol(transformedData))
  var_explained <- data.frame(PCs = 1:length(dat_svd$d), var = cumsum(dat_svd$d/sum(dat_svd$d)))
  print(var_explained)
  n_pcs <- predict(lm(PCs~poly(var,2,raw=T),var_explained[1:nrow(var_explained),]),
                 newdata = data.frame(var=1))
  dat_reduced <- data.table(dat_svd$u[,n_pcs] %*% diag(dat_svd$d[n_pcs], length(n_pcs), length(n_pcs)) %*% t(dat_svd$v[,n_pcs]))
  return(list(var_explained, data.table(dt[,list(ID_REF,IDENTIFIER)],
                                        dat_reduced, 
              key = c("ID_REF","IDENTIFIER"))))
}

result <- preprocessData(res[[2]])
gc()

xyplot(var~PCs,result$S.cumexplvar,panel = function(x,y,...){
  panel.xyplot(x,y,...)
  nfit <- round(log2(length(result$S.cumexplvar$PCs)))
  tmp <- predict(lm(PCs~poly(var,2,raw=T),result$S.cumexplvar[1:nfit,]),data.frame(var=y))
  res <- predict(lm(PCs~poly(var,2,raw=T),result$S.cumexplvar[1:nfit,]),data.frame(var=1))
  panel.xyplot(tmp,y, type = "l")
  panel.abline(v = res)
})


res[,lapply(.SD, function(x) is.na(x))), .SDcols = names(res)]


colsums <- res[,lapply(.SD,function(x) sum(is.na(x))),.SDcols = names(res)]
rowsums <- res[,lapply(.SD, is.na), .SDcols  = names(res)][,rowSums(.SD),by=ID_REF, .SDcols = grep("GSM",names(res))]

res[IDENTIFIER=="--Control"]

scale2 <- scale(t(res[,lapply(.SD,function(x)scale(log2(x))),.SDcols=grep("GSM",names(res))]))
resscaled <- res[,lapply(.SD,function(x)scale(log2(x))),.SDcols=grep("GSM",names(res))]
res[complete.cases(res)][,lapply(.SD,mean),.SDcols = grep("GSM",names(res))]
res[complete.cases(res)][,.SDcols = grep("GSM",names(res))]


dim(res[complete.cases(res)])
dim(res)


res[,.SDcols = grep("GSM",names(res))]

print(res[complete.cases(res)][!complete.cases(res)],nrow=100)


res[[2]][,lapply(.SD,mean), .SDcols = names(res[[2]])]
test <- res[[2]][complete.cases(res[[2]])]
