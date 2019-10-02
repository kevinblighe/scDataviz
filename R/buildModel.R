buildModel <- function(
  data,
  group = NULL,
  groupfactors = unique(data$metadata[,group]),
  cluster1 = NULL,
  cluster2 = NULL,
  markers = colnames(data$expression))
{
  modeling <- data$expression[,which(colnames(data$expression) %in% markers)]

  if (!is.null(group1)) {
    modeling$group <- factor(data$metadata[,group], levels = groupfactors)

    # build a predictive model of OSSTAT
    f <- as.formula(
      paste0('group ~ ',
        paste0(markers, collapse = "+")))

    # use step-wise regression with AIC in order to reduce model in size and
    # select 'best' predictors
    require(arm)
    modeling[is.na(modeling)] <- 0
    null <- bayesglm(group ~ 1, data = modeling, family = binomial(link = 'logit'))
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