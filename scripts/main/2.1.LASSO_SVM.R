library(data.table)
library(dplyr)
library(pROC)
library(caret)
set.seed(456)
library(PreProcess)

directory = "/path/to/data_dir"

input_dir=paste0(directory, "/data/dataInput")
output_dir=paste0(directory, "/data/test/Model/LASSO")
if (!dir.exists(output_dir)) {dir.create(output_dir, recursive=T)}


dt=as.data.frame(fread(paste0(input_dir, "/bigTable.mergedAll.tsv")))
row.names(dt)=dt$probeID
# tune data
sampleSheet=read.table(paste0(directory, "/config/train.tsv"),sep="\t", header=T)
# cpgList=read.table(paste0(directory, "/data/train/DMP/TB-HC/DMPs.hg38.tsv"), sep="\t", header=T)

SigCpG=read.table(paste0(directory, "/data/train/DMP/TB-HC/LASSO.txt"), sep="\t", header=T)$probeID

# cpgList=cpgList[cpgList$probeID%in%SigCpG, ]

df=dt[dt$probeID%in%SigCpG, sampleSheet$Accession]

tdf=data.frame(t(df))
tdf$sample_ID=rownames(tdf)
tdf$Class=merge(tdf, sampleSheet, by.x="sample_ID", by.y="Accession")$Class
tuneDT=arrange(tdf, desc(Class))
tuneDT=tuneDT[,-grep("sample_ID", names(tuneDT))]
# tuneDT$Class=ifelse(tuneDT$Class=="HC", 0, 1)

# validation data
ValSample=read.table(paste0(directory, "/config/test.tsv"),sep="\t", header=T)
df=dt[dt$probeID%in%SigCpG, ValSample$Accession]
tdf=data.frame(t(df))
tdf$sample_ID=rownames(tdf)
tdf$Class=merge(tdf, ValSample, by.x="sample_ID", by.y="Accession")$Class
valDT=arrange(tdf, desc(Class))
valDT=valDT[,-grep("sample_ID", names(valDT))]
# valDT$Class=ifelse(valDT$Class=="HC", 0, 1)

# scale <- preProcess(tuneDT[, 1:5], method = "scale")
# valDT <- predict(scale, valDT)
# tuneDT<-predict(scale, tuneDT)

### Model development
# Run algorithms using 10-fold cross validation
trainControl=trainControl(method="repeatedcv", number=10, repeats=3)

## the tuning Grid
rfGrid <- expand.grid(C = 1:(ncol(tuneDT)-1))
seeds.svm <- vector(mode="list", length=31)

## the number of k in knn 10, the same number of C in rf
for(i in 1:30) seeds.svm[[i]] <- sample.int(n=1000, 25)
# for the last model
seeds.svm[[31]] <- sample.int(1000, 1)

# no sampling for imbalance class
print("the trControl for caret 10 cross-fold 3 repeats no smote")
fitControl.svm <- trainControl(
  method="repeatedcv",
  number=10,
  repeats=3,
  allowParallel=TRUE,
  classProbs=TRUE,
  seeds=seeds.svm,
  summaryFunction=twoClassSummary,
  verboseIter=FALSE)

print("no smote sampling, Model1 rf")
rf.m1 <- train(Class ~ ., 
               data=tuneDT,
               method="svmLinear",
               metric="ROC",
               trControl=fitControl.svm,
               tuneGrid=rfGrid)

pdf(file.path(output_dir, "/SVM_Cost_CrossValidation.pdf"), width=3, height=3)
plot(rf.m1)
dev.off()

pre.m1s <- predict(rf.m1, tuneDT, type="prob")
# confusion matrix for training
pre.m1s$pred <- factor(ifelse(pre.m1s$HC >= .5, "HC", "TB"))
cm <- confusionMatrix(pre.m1s$pred, as.factor(tuneDT$Class))
m1_training_all <- as.data.frame(cm$byClass)
colnames(m1_training_all) <- "nosmote_training"
# roc for training
roc.m1s <- roc(tuneDT$Class, pre.m1s$HC)
pROC::auc(roc.m1s)
pROC::ci.auc(roc.m1s)
coords(roc.m1s, "best", ret="all")

pdf(file.path(output_dir, "/SVM_AUC_testDT.pdf"))
ciobj <- ci.se(roc.m1s, specificities=seq(0, 1, l=length(roc.m1s$sensitivities)))
dat.ci <- data.frame(x = as.numeric(rownames(ciobj)),
                     lower = ciobj[, 1],
                     upper = ciobj[, 3])
auc <- paste0("AUC = " ,round(roc.m1s$auc, 2))

### Plot
ggroc(roc.m1s, size = 2, color = "darkblue", legacy.axes = TRUE) + 
  coord_equal() + 
  theme_minimal() + 
  geom_abline(slope = 1, intercept = 0, size =1,
              linetype = "dashed", alpha=0.7, color = "grey") + 
  geom_ribbon(data = dat.ci, aes(x = 1-x, ymin = lower, ymax = upper), 
              fill = "steelblue", alpha= 0.2) + 
  labs(x = "1 - Specificity", y = "Sensitivity") +
  annotate("text", x = 1, y = 0, label = paste0(auc), size = 10,
           colour = "black" , vjust = "inward", hjust = "inward") + 
  theme(text = element_text(size = 15)) 

dev.off()

you.m1s <- coords(roc.m1s, x="best", best.method="youden", ret="all")
m1_train <- as.data.frame(you.m1s)
m1_train["AUC",] <- roc.m1s$auc
m1_train["mtry",] <- rf.m1$best
colnames(m1_train) <- "nosmote_training"

pre.m1s <- predict(rf.m1, valDT, type="prob")
# confusion matrix for testing
pre.m1s$pred <- factor(ifelse(pre.m1s$HC >= .5, "HC", "TB"))
cm <- confusionMatrix(pre.m1s$pred, as.factor(valDT$Class))
m1_testing_all <- as.data.frame(cm$byClass)
colnames(m1_testing_all) <- "nosmote_testing"
# roc for testing
roc.m1s <- roc(valDT$Class, pre.m1s$HC)
you.m1s <- coords(roc.m1s, x="best", best.method="youden", ret="all")
m1_test <- as.data.frame(you.m1s)
m1_test["AUC",] <- roc.m1s$auc
m1_test["C",] <- rf.m1$best
colnames(m1_test) <- "nosmote_testing"
# plot ROC for testing
rf.m1.ROC <- roc(predictor=pre.m1s$HC, response=valDT$Class)
pROC::auc(rf.m1.ROC)
pROC::ci.auc(rf.m1.ROC)
coords(rf.m1.ROC, "best", ret="all")

auc <- pROC::auc(rf.m1.ROC)
auc.ci <- pROC::ci.auc(rf.m1.ROC)
pdf(file.path(output_dir, "/SVM_AUC_valDT.pdf"))
ciobj <- ci.se(rf.m1.ROC, specificities=seq(0, 1, l=length(roc.m1s$sensitivities)))
dat.ci <- data.frame(x = as.numeric(rownames(ciobj)),
                     lower = ciobj[, 1],
                     upper = ciobj[, 3])
auc <- paste0("AUC = " ,round(rf.m1.ROC$auc, 2))

### Plot
ggroc(rf.m1.ROC, size = 2, color = "darkblue", legacy.axes = TRUE) + 
  coord_equal() + 
  theme_minimal() + 
  geom_abline(slope = 1, intercept = 0, size =1,
              linetype = "dashed", alpha=0.7, color = "grey") + 
  geom_ribbon(data = dat.ci, aes(x = 1-x, ymin = lower, ymax = upper), 
              fill = "steelblue", alpha= 0.2) + 
  labs(x = "1 - Specificity", y = "Sensitivity") +
  annotate("text", x = 1, y = 0, label = paste0(auc), size = 10,
           colour = "black" , vjust = "inward", hjust = "inward") + 
  theme(text = element_text(size = 15)) 
dev.off() 
