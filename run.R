library(GEOmetadb)
# library(GEOquery)
library(plyr)
library(data.table)
library(caret)
library(doMC)
library(exprt)
registerDoMC(18)

# getSQLiteFile(destdir ="../genePhenotypePrediction/" )
con <- dbConnect(SQLite(),'../genePhenotypePrediction/GEOmetadb.sqlite')
gdsquery <- dbGetQuery(con,'SELECT sample_organism, gds, sample_count, gpl, value_type from gds WHERE (sample_organism="Homo sapiens" OR sample_organism="Mus musculus") AND (value_type="count" OR value_type="transformed count")')
dbDisconnect(con)

#
# GAD <- fread("/scratch/jfalck/usr/rprojects/exprt/data/GADCDC/GADCDC_data.tsv",
#              sep = "\t", header = T, )
#


gdsquery <- data.frame(lapply(gdsquery,function(x) gsub(" ","_",x)))
gdsquery$sample_count <- as.integer(gdsquery$sample_count)
gds <- Filter(function(x){sum(x$sample_count) >= 50},split(gdsquery, list(gdsquery$sample_organism, gdsquery$gpl, gdsquery$value_type),))


initializeExpr <- function(gds = "list"){
  lapply(gds, function(o){
    new("Expr", organism = gsub(" ","_",as.character(unique(o$sample_organism))),
        valueType = gsub(" ","_",as.character(unique(o$value_type))),
        gpl = as.character(unique(o$gpl)),
        gds = as.character(unique(o$gds)),
        sampleCounts = as.integer(o$sample_count))
  })
}



# cleft_lip_genes <- toupper(GAD[DISEASE %like% "cleft lip",unique(GENE)])



test <- initializeExpr(gds)


analysis <- function(i){
  plot(i$var_explained,
       ylim = c(0,1))
  abline(v = i$n_pcs)
}
#
# library(pROC)
#
#
# ho <- sample(cleft_lip_genes,80,replace = F)
# train <- cleft_lip_genes[!cleft_lip_genes %in% ho]
#
# testl <- c("Homo_sapiens.GPL97.transformed_count","Homo_sapiens.GPL2895.count",
#            "Homo_sapiens.GPL571.count","Homo_sapiens.GPL571.transformed_count",
#            "Homo_sapiens.GPL96.transformed_count","Homo_sapiens.GPL97.count")

tbd <- file.exists(paste(gsub("[.]","_",names(test)),"rds",sep="."))
test <- test[!tbd]

e <- new.env()
e$dat <- ""
doforall <- function(ds){
  test[[ds]]$load()
  test[[ds]]$merge(fbind = F)
  # test[[ds]]$impute()
  test[[ds]]$reduce()
  test[[ds]]$exprSvd$save()
#   e$dat <<- readRDS(paste(names(test)[ds],"rds.gz",sep = "."),compress="gzip")
#   e$dat$reduce(.6)
#   # test[[ds]]$save(names(test)[ds])
#   e$dat$pred(train,666)
#   e$dat$prediction[,.(IDENTIFIER,
#                           "prediction" = prediction,
#                           "pcs" = e$dat$n_pcs)][]
}

res <- l_ply(seq_along(test)[95:length(test)],doforall)

#
# rocfun <- function(i){
#   plot(roc(res[[i]]$IDENTIFIER %in% ho, res[[i]]$prediction))
# }
#
# total <- 20
# for(i in 1:total){
#   Sys.sleep(0.1)
#   cat(i)
#   # update GUI console
#   flush.console()
# }
#
#
#
#
# tdt <- data.table::rbindlist(res)
# tdt <- tdt[,.("prediction" = mean(prediction*pcs)/sum(pcs)),by=IDENTIFIER]
#
# plot(roc(tdt$IDENTIFIER %in% ho, tdt[,(prediction)]))



mdat <- test[[4]]$gdsMissing[,3:ncol(test[[4]]$gdsMissing),with=F]
