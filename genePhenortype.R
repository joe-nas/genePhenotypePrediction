library(GEOmetadb)
library(GEOquery)
library(plyr)
library(data.table)
library(caret)
# download database
# if(!file.exists('GEOmetadb.sqlite')) getSQLiteFile()

# library(doParallel)
# nodes <- detectCores()
# cl <- makeCluster(6)
# registerDoParallel(cl)
#stopCluster(cl)
library(doMC)
registerDoMC(10)



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
# library(doMC)
# registerDoMC(4)
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
  return(dt)
}

softpath <- function(x,i,workingdir,gdsline){
  x$gds[[i]] <- sprintf("%s/%s/%s/%s/%s.soft.gz", workingdir, 
                        gdsline$sample_organism, gdsline$gpl, 
                        gdsline$value_type, gdsline$gds)
  invisible()
}

gpldtpath <- function(workingdir, gdsline){
  sprintf("%s/%s/%s_%s_%s_table.csv", workingdir, gdsline$sample_organism, 
          gdsline$sample_organism, gdsline$gpl, gdsline$value_type)
}

savedt <- function(gpldtpath, dt){
  write.table(dt ,file = gpldtpath, col.names = T, row.names = F, sep = ",")
}


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

l_ply(seq_along(e$gds)[1:10], function(i){
  modify(e,i,llply(e$gds[[i]], soft2dt, .parallel = T))
  invisible()
})


l_ply(seq_along(e$gds)[1:10], function(i){
  modify(e,i, Reduce(function(k,j) k[j], Map(function(z) return(z) ,e$gds[[i]])))
  invisible()
},.parallel = F)

modifyl(e,"gds", Filter(function(x) dim(x)[2]>50, e$gds))



dtProcessing <- function(e,i){
  # imputing missing values
  preproc <- preProcess(e$gds[[i]][,3:ncol(e$gds[[i]]),with=F], 
                        method = c("center","scale","medianImpute"))
  e$gds[[i]][,3:ncol(e$gds[[i]]) := data.table(predict(preproc,e$gds[[i]][,3:ncol(e$gds[[i]]),with=F])), with=F]
  
  #extracting right singular values
  dt_svd <- svd(e$gds[[i]][,3:ncol(e$gds[[i]]), with = F], nv = ncol(e$gds[[i]])-2, nu = 0)
  var_explained <- data.table(PCs = 1:length(dt_svd$d), Var = cumsum(dt_svd$d/sum(dt_svd$d)))
  slope <- diff(var_explained[,Var])/diff(var_explained[,PCs])
  n_pcs <- max(which(slope > quantile(slope,.90)))
  e$gds[[i]][e$gds[[i]][,1:2, with=F],data.table((as.matrix(e$gds[[i]][,3:ncol(e$gds[[i]]),with=F]) %*% t(dt_svd$v))[,1:n_pcs])]
  write.table(e$gds[[i]],file = names(e$gds[i]), sep = ",",col.names = F,row.names = T,quote = F)
  return(list(var_explained))
}





testres<-lapply(seq_along(e$gds)[1:4],function(i) dtProcessing(e,i))
llply(seq_along(e$gds),function(i)names(e$gds[i]))


dat <- testres[[4]][[1]]
diffy <- diff(dat[,Var])
slope <- diff(dat[,Var])/diff(dat[,PCs])

plot(dat)
abline(v=max(which(slope > quantile(slope,.90))))



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

