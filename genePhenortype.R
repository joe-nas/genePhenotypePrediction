library(GEOmetadb)
library(GEOquery)
library(plyr)
library(data.table)
library(caret)
library(doMC)
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
expression_data <- llply(rev(seq_along(gds)), function(i){
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
              verbose = F)[, .SD,, c("ID_REF","IDENTIFIER")]
  cat(".")
  return(dt)
}

softpath <- function(x,i,workingdir,gdsline){
  x$gds[[i]] <- sprintf("%s/%s/%s/%s/%s.soft.gz", workingdir, gdsline$sample_organism, gdsline$gpl, gdsline$value_type, gdsline$gds)
  invisible()
}

gpldtpath <- function(workingdir, gdsline){
  sprintf("%s/%s/%s_%s_%s_table.csv", workingdir, gdsline$sample_organism, 
          gdsline$sample_organism, gdsline$gpl, gdsline$value_type)
}

savedt <- function(gpldtpath, dt){
  write.table(dt ,file = gpldtpath, col.names = T, row.names = F, sep = ",")
}

<<<<<<< HEAD
registerDoMC(3)



e <- new.env()
e$gds <- gds

modify <- function(x, i, arg){
  x$gds[[i]] <- arg
  invisible()
}

modifyl <- function(x,list_name, arg){
  x[[list_name]] <- arg
  invisible()
}

l_ply(seq_along(e$gds), function(i){
  softpath(e,i,".",e$gds[[i]])
  modify(e, i, e$gds[[i]][file.exists(e$gds[[i]])])
  invisible()
  })
length(e$gds)
=======
registerDoMC(1)
Sys.setenv(OPENBLAS_NUM_THREADS = 4)
Sys.getenv("OPENBLAS_NUM_THREADS")
preprocessData <- function(dt){
  require(caret)
  dt_dim <- dim(dt)
  PreProcValues <- preProcess(dt[,(.SD),.SDcols = grep("GSM",names(dt))], 
                              method = c("center","scale","medianImpute"))
  svd_data <- svd(
    predict(PreProcValues, dt[,(.SD),.SDcols = grep("GSM",names(dt))]),
    nu = dt_dim[1], nv = dt_dim[2])
  var_explained <- data.frame(PCs = 1:length(svd_data$d), Var = cumsum(svd_data$d/sum(svd_data$d)))
  nfit <- round(log2(length(var_explained$PCs)))
  so_model <- lm(PCs~poly(Var,2,raw=T), data = var_explained[1:nfit,])
  n_pcs <- predict(so_model, newdata = data.frame(Var=1))
  n_pcs <- 1:round(n_pcs)
  dat_reduced <- data.table(svd_data$u[,n_pcs] %*% diag(svd_data$d[n_pcs], length(n_pcs), length(n_pcs)))
  return(list(var_explained = var_explained, 
              so_model = so_model,
              data.table(
                dt[,list(ID_REF,IDENTIFIER)], 
                dat_reduced, 
                key = c("ID_REF","IDENTIFIER")) 
              ))
}


library(lattice)
xyplot(Var~PCs, data = result$var_explained,
       so_model = result$so_model,
       panel = function(x,y,so_model,...){
         target_pc <- round(predict(so_model,data.frame(Var=1)))
         xso <- predict(so_model, newdata = data.frame(Var = result$var_explained$Var))
         panel.xyplot(x,y,...)
         panel.xyplot(xso,y, ..., type = "l")
         panel.abline(v = target_pc)
})






e <- new.env()
e$gds <- gds

modify <- function(x, i, arg){
  x$gds[[i]] <- arg
  invisible()
}

modifyl <- function(x,list_name, arg){
  x[[list_name]] <- arg
  invisible()
}

l_ply(seq_along(e$gds), function(i){
  softpath(e,i,".",e$gds[[i]])
  modify(e, i, e$gds[[i]][file.exists(e$gds[[i]])])
  invisible()
  })
length(e$gds)

modifyl(e,"gds",Filter(length, e$gds))

registerDoMC(12)
l_ply(seq_along(e$gds)[1:5], function(i){
  modify(e,i,llply(e$gds[[i]], soft2dt, .parallel = T))
  invisible()
})


l_ply(seq_along(e$gds)[1:5], function(i){
  modify(e,i, Reduce(function(x,y) x[y], e$gds[[i]]))
  invisible()
})

modifyl(e,"gds", Filter(function(x) dim(x)[2]>50, e$gds))
parent.env(e)



library(caret)
testdat <- e$gds[[1]]

Sys.setenv(OMP_NUM_THREADS=6,OPENBLAS_NUM_THREADS=6)

dtProcessing <- function(dt){
  preproc <- preProcess(dt[,3:ncol(dt),with=F], 
                        method = c("center","scale","medianImpute"))
  dt[,3:ncol(dt) := data.table(predict(preproc,dt[,3:ncol(dt),with=F])), with=F]
  dt_svd <- svd(as.matrix(dt[,3:ncol(dt), with = F]),nu = nrow(dt),nv = ncol(dt)-2)
  var_explained <- data.table(PCs = 1:length(dt_svd$d), Var = cumsum(dt_svd$d/sum(dt_svd$d)))
  pcs_fit <- round(log2(length(var_explained$PCs)))
  fit <- lm(PCs~poly(Var,2,raw=T), data = var_explained[1:pcs_fit,])
  n_pcs <- predict(fit, newdata = data.frame(Var=1))
  n_pcs <- 1:round(n_pcs)
  dt[,3:ncol(dt) := NULL]
  dt[data.table(dt_svd$u[,n_pcs] %*% diag(dt_svd$d[n_pcs], length(n_pcs), length(n_pcs)))]
  
}

dtProcessing(testdat)


testpreproc <- preProcess(testdat[,3:ncol(testdat),with=F], 
                          method = c("center","scale","medianImpute"))
testmat <- as.matrix(testdat[,3:ncol(testdat) := data.table(predict(testpreproc,testdat[,3:ncol(testdat),with=F])), with=F][,3:ncol(testdat),with=F])

mysvd <- svd(testmat, nu = nrow(testmat), nv = ncol(testmat))












grep("GSM",names(testdat),value = T)

l_ply(seq_along(e$gds)[1:3], function(i){
  modify(e,i,preprocessData(e$gds[[i]]))
  invisible()
}, .progress = "text")










>>>>>>> 065ca5314f42663fa81bc4c37b007ffd6c0096a6

modifyl(e,"gds",Filter(length, e$gds))

# registerDoMC(12)
l_ply(seq_along(e$gds)[1:5], function(i){
  modify(e,i,llply(e$gds[[i]], soft2dt, .parallel = T))
  invisible()
})


l_ply(seq_along(e$gds)[1:5], function(i){
  modify(e,i, Reduce(function(x,y) x[y], e$gds[[i]]))
  invisible()
})

modifyl(e,"gds", Filter(function(x) dim(x)[2]>50, e$gds))
parent.env(e)



library(caret)
testdat <- e$gds[[3]]

# Sys.setenv(OMP_NUM_THREADS=6,OPENBLAS_NUM_THREADS=6)

dtProcessing <- function(dt){
  # imputing missing values
  preproc <- preProcess(dt[,3:ncol(dt),with=F], 
                        method = c("center","scale","medianImpute"))
  dt[,3:ncol(dt) := data.table(predict(preproc,dt[,3:ncol(dt),with=F])), with=F]
  
  #extracting right singular values
  dt_svd <- svd(dt[,3:ncol(dt), with = F], nv = ncol(dt)-2, nu = 0)
  var_explained <- data.table(PCs = 1:length(dt_svd$d), Var = cumsum(dt_svd$d/sum(dt_svd$d)))
  
  
  pcs_fit <- round(log2(length(var_explained$PCs)))
  fit <- lm(PCs~poly(Var,2,raw=T), data = var_explained[1:pcs_fit,])
  n_pcs <- predict(fit, newdata = data.frame(Var=1))
  n_pcs <- n_pcs
  print(n_pcs)
  #dt[dt[,1:2, with=F],data.table((as.matrix(dt[,3:ncol(dt),with=F]) %*% t(dt_svd$v))[,1:n_pcs])]
  return(list(var_explained,fit))
  
}

testres<-dtProcessing(testdat)


library(lattice)
xyplot(Var~PCs, data = result$var_explained,
       so_model = result$so_model,
       panel = function(x,y,so_model,...){
         target_pc <- round(predict(so_model,data.frame(Var=1)))
         xso <- predict(so_model, newdata = data.frame(Var = result$var_explained$Var))
         panel.xyplot(x,y,...)
         panel.xyplot(xso,y, ..., type = "l")
         panel.abline(v = target_pc)
       })

