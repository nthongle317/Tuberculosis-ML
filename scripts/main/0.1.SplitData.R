library(caret)
library(data.table)

set.seed(1)

dt=read.table("inputSample.csv", sep=",", header=T)

### setup train and test data
index=createDataPartition(dt$Class, p=0.70, list=FALSE,  times = 1)
# select 30% of the data for testing
testDT=dt[-index,]
# use the remaining 70% of data to training and testing the models
trainDT=dt[index,]

table(trainDT$Class)
table(testDT$Class)

write.table(trainDT, "train.tsv", sep="\t", col.names = T, quote=F, row.name=F)
write.table(testDT, "test.tsv",  sep="\t", col.names = T, quote=F, row.name=F)