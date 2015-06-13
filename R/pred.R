#' @export
Expr$methods(
  pred = function(gene_list){
    y <- as.factor(gdsRed[,toupper(IDENTIFIER)] %in% toupper(gene_list))
    levels(y) <- gsub("TRUE", "yes", levels(y))
    levels(y) <- gsub("FALSE", "no", levels(y))
    x <- quote(gdsRed[,c(3:ncol(gdsRed)),with = F])
    ts <- downSample(eval(x),y)

    tc <- trainControl(classProbs = T, method = "cv",
                       number = 10, allowParallel = T)

    fit <- train(Class ~ ., data = ts,
                  method = "rf", trControl = tc)
    prediction <<- predict(fit, newdata =  eval(x),
                           type = "prob")
    prediction <<- data.frame(IDENTIFIER = gdsRed$IDENTIFIER,
                              prediction = prediction)
  }
)
