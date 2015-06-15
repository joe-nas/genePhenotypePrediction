#' @export
Expr$methods(
  pred = function(gene_list,seed = NA){
    y <- as.factor(gdsRed[,toupper(IDENTIFIER)] %in% toupper(gene_list))
    levels(y) <- gsub("TRUE", "yes", levels(y))
    levels(y) <- gsub("FALSE", "no", levels(y))
    x <- quote(gdsRed[,c(3:ncol(gdsRed)),with = F])
    set.seed(666)
    ts <- downSample(eval(x),y)


    tc <- trainControl(classProbs = T, method = "cv",
                       number = 10, allowParallel = T)

    cat(sprintf("fitting\n"))
    set.seed(seed)
    fit <<- train(Class ~ ., data = ts,
                  method = "rf", trControl = tc)
    cat(sprintf("predicting\n"))
    prediction <<- data.table(
      gdsData[,.(IDENTIFIER)],
      predict(fit, newdata = eval(x), type = "prob")
      )[,.("prediction" = mean(yes)), by=IDENTIFIER]


    invisible(prediction)
  }
)
