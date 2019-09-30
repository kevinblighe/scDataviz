buildModel <- function()
{
  modeling <- read.csv("modeling.txt", sep = '\t', stringsAsFactors = FALSE)

  # removes markers that were used for gating (?)
  modeling <- modeling[,-which(colnames(modeling) %in% c('CD45','CD4','CD3'))]

  # retrospectively remove markers of high std error
  modeling <- modeling[,-which(colnames(modeling) %in% c('IFN.g','C3aNeo_i'))]

  modeling <- data.frame(
    Sample = modeling$Sample,
    asinh(
      data.frame(
        data.matrix(modeling[,2:ncol(modeling)])/5)))

  modeling$Sample <- factor(modeling$Sample, levels = c('HD', 'Patient'))

  # determine markers 'co-expressed' with 3 markers of interest
  mat <- modeling
  scale01 <- function(x){(x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))}
  mat[,2:ncol(mat)] <- scale01(mat[,2:ncol(mat)])
  mat[is.na(mat)] <- 0
  C5AR1_top20 <- order(mat[which(mat$Sample == 'Patient'),'C5aR1_i'], decreasing = TRUE)[1:10]
  C3C3B_top20 <- order(mat[which(mat$Sample == 'Patient'),'C3.C3b_i'], decreasing = TRUE)[1:10]
  C5C5B_top20 <- order(mat[which(mat$Sample == 'Patient'),'C5.C5b_i'], decreasing = TRUE)[1:10]
  topten <- names(sort(apply(mat[as.numeric(names(which(table(c(C5AR1_top20, C3C3B_top20, C5C5B_top20)) == 2))),2:ncol(mat)], 2, mean), decreasing = TRUE)[1:10])
  mat <- data.matrix(mat[,topten])
  rownames(mat) <- modeling$Sample
  heatmap(mat)

  # build a predictive model of OSSTAT
  f <- as.formula(
    paste0("Sample ~ ",
      paste0(colnames(modeling[,2:ncol(modeling)]), collapse = "+")))

  #f_topten <- as.formula(paste0("Sample ~ ", topten))
  #f <- f_topten

  # use step-wise regression with AIC in order to reduce model in size and
  # select 'best' predictors
  require(arm)
  modeling[is.na(modeling)] <- 0
  null <- bayesglm(Sample ~ 1, data = modeling, family = binomial(link = 'logit'))
  full <- bayesglm(f, data = modeling, family = binomial(link = 'logit'))

  require(MASS)
  forward <- stepAIC(null, scope=list(lower=null, upper=full),
    direction="forward", k = log(nrow(modeling)))
  backward <- stepAIC(full,
    direction="backward", k = log(nrow(modeling)))
  both <- stepAIC(null, scope=list(lower=null, upper=full),
    direction="both", k = log(nrow(modeling)))

  # remove markers with unusual error and/or odds ratio CIs
  anova(full, test="Chisq")
  anova(both, test="Chisq")
  #both <- update(both, . ~ . - monopoiesis_HLA)
  anova(forward, test="Chisq")
  #forward <- update(forward, . ~ .- monopoiesis_HLA)
  anova(backward, test="Chisq")
  #backward <- update(backward, . ~ . - events_45_ - events_gran - events_45_1 - monopoiesis_HLA - monopoiesis_HLA_Lympho - myeloid_progenitors_CD117 - erythropoiesis_MFI_CD105_CD105_Ery - erythropoiesis_CD235a)

  # write ouy tables
  a <- data.frame(rownames(coefficients(summary(full))), coefficients(summary(full)), exp(cbind("Odds ratio" = coef(full), confint.default(full, level = 0.90))))[-1,]
  colnames(a) <- c("Marker", "Beta", "Std. Error", "Z value", "P value", "OR", "Lower CI", "Upper CI")
  a$Marker <- gsub("_", " ", a$Marker)
  write.table(a, "FullModel.tsv", sep = "\t", row.names = FALSE, quote = FALSE, dec = ".")
  a <- data.frame(rownames(coefficients(summary(both))), coefficients(summary(both)), exp(cbind("Odds ratio" = coef(both), confint.default(both, level = 0.90))))[-1,]
  colnames(a) <- c("Marker", "Beta", "Std. Error", "Z value", "P value", "OR", "Lower CI", "Upper CI")
  a$Marker <- gsub("_", " ", a$Marker)
  write.table(a, "Model1.tsv", sep = "\t", row.names = FALSE, quote = FALSE, dec = ".")
  a <- data.frame(rownames(coefficients(summary(forward))), coefficients(summary(forward)), exp(cbind("Odds ratio" = coef(forward), confint.default(forward, level = 0.90))))[-1,]
  colnames(a) <- c("Marker", "Beta", "Std. Error", "Z value", "P value", "OR", "Lower CI", "Upper CI")
  a$Marker <- gsub("_", " ", a$Marker)
  write.table(a, "Model2.tsv", sep = "\t", row.names = FALSE, quote = FALSE, dec = ".")
  a <- data.frame(rownames(coefficients(summary(backward))), coefficients(summary(backward)), exp(cbind("Odds ratio" = coef(backward), confint.default(backward, level = 0.90))))[-1,]
  colnames(a) <- c("Marker", "Beta", "Std. Error", "Z value", "P value", "OR", "Lower CI", "Upper CI")
  a$Marker <- gsub("_", " ", a$Marker)
  write.table(a, "Model3.tsv", sep = "\t", row.names = FALSE, quote = FALSE, dec = ".")

  # perform 10-fold cross-validation on final models to test their utility
  require(boot)
  cv.glm(modeling, full, K=10)$delta
  cv.glm(modeling, both, K=10)$delta
  cv.glm(modeling, forward, K=10)$delta
  cv.glm(modeling, backward, K=10)$delta

  # perform ROC analysis and generate a single ROC curve for all models
  pdf("ROC.pdf", width=8, height=8)
    par(mar = c(5,5,5,5))

    require(pROC)
    roc <- roc(modeling$Sample, fitted(full), smooth = FALSE)
    plot.roc(roc,
      grid = TRUE,
      grid.lwd = 2,
      col = "black",
      main = "ROC analysis - unstimulated data\nReduced models selected via AIC criterion, and cross validated 10-fold",
      sub = "")
    text(0.255, 0.35,
      paste("Full model [AUC (L95; U95)]:\n",
        paste0(round(ci.auc(roc)[2], 3), " (", round(ci.auc(roc)[1], 3), "; ", round(ci.auc(roc)[3], 3), ")"), sep=""),
      cex=1.0, col = "black")

    roc <- roc(modeling$Sample, fitted(both), smooth = FALSE)
    lines.roc(roc, grid=TRUE, grid.lwd=2, col = "red2")
    text(0.255, 0.25,
      paste("Reduced model 1 [AUC (L95; U95)]:\n",
        paste0(round(ci.auc(roc)[2], 3), " (", round(ci.auc(roc)[1], 3), "; ", round(ci.auc(roc)[3], 3), ")"), sep=""),
      cex=1.0, col = "red2")

    roc <- roc(modeling$Sample, fitted(forward), smooth = FALSE)
    lines.roc(roc, grid=TRUE, grid.lwd=2, col = "royalblue")
    text(0.255, 0.15,
      paste("Reduced model 2 [AUC (L95; U95)]:\n",
        paste0(round(ci.auc(roc)[2], 3), " (", round(ci.auc(roc)[1], 3), "; ", round(ci.auc(roc)[3], 3), ")"), sep=""),
      cex=1.0, col = "royalblue")

    roc <- roc(modeling$Sample, fitted(backward), smooth = FALSE)
    lines.roc(roc, grid=TRUE, grid.lwd=2, col = "forestgreen")
    text(0.255, 0.05,
      paste("Reduced model 3 [AUC (L95; U95)]:\n",
        paste0(round(ci.auc(roc)[2], 3), " (", round(ci.auc(roc)[1], 3), "; ", round(ci.auc(roc)[3], 3), ")"), sep=""),
      cex=1.0, col = "forestgreen")
  dev.off()
}