library(tidyverse)
library(stringr)
library(dplyr)
library(ggplot2)
library(factoextra)
library(NbClust)
library(caret)
library(caTools)
library(corrplot)
library(cluster)
library(e1071)
library(lattice)
##### read files 
bcancer_raw <- read.csv("~/Desktop/DANA 4840/project/77_cancer_proteomes_CPTAC_itraq.csv")
patient <- read.csv("~/Desktop/DANA 4840/project/clinical_data_breast_cancer.csv")
pam50 <- read.csv("~/Desktop/DANA 4840/project/PAM50_proteins.csv")


############### Combining Datasets ##################
bcancer <- data.frame(t(bcancer_raw[-1]))
colnames(bcancer) <- bcancer_raw[,1]
bcancer <- bcancer[3:(nrow(bcancer)-3), ]
TCGA <- rownames(bcancer)
bcancer$TCGA <- TCGA
bcancer <- tibble::rowid_to_column(bcancer, "ID")

#write.csv(bcancer,"~/Desktop/DANA 4840/project/bcancer.csv")

# cleaning for TCGA in two csv file

bcancer$TCGA<-str_sub(bcancer$TCGA, end=-8)
bcancer$TCGA<-gsub('[[:punct:] ]+','',bcancer$TCGA)
patient$Complete.TCGA.ID <- str_sub(patient$Complete.TCGA.ID, start = 5)
patient$Complete.TCGA.ID <- gsub('[[:punct:] ]+','',patient$Complete.TCGA.ID)
colnames(patient)[1]<-'TCGA'

# Drop columns has more than 50% NAs 
bcancer <- bcancer[, which(colMeans(!is.na(bcancer)) > 0.5)] 
#mice_data<-mice(bcancer,m=1,maxit=2,seed=333)
#bcancer<-complete(mice_data,1)


#install.packages('varhandle')
library('varhandle')
bcancer[,1:11512]<- unfactor(bcancer[,1:11512])
for(i in 2:11512){
  bcancer[is.na(bcancer[,i]), i] <- median(bcancer[,i], na.rm = TRUE)
}
View(bcancer)

# Left-join two dataset
combined_bcancer <- left_join(bcancer, patient, by = 'TCGA', copy = FALSE)

View(combined_bcancer)
write.csv(combined_bcancer,"~/Desktop/DANA 4840/project/combined_bcancer.csv")




######################### Start From Here #####################


df <- read.csv('../Group Project/combined_bcancer (1).csv')
PAM50 <- read.csv('../Group Project/PAM50_proteins.csv')

df <- df[which(df$Gender != 'MALE'), ] #Removing Males

df$HER2.Final.Status <- ifelse(df$HER2.Final.Status == 'Equivocal'|df$HER2.Final.Status == 'Positive', 'Positive', 'Negative')
### Matching ids with PAM50 dataset

id_PAM50 <- PAM50$RefSeqProteinID

id_bcancer <- colnames(df)

new_df_feature_important <- df[which(id_bcancer %in% id_PAM50)] %>% scale()

dim(new_df_feature_important)

head(new_df_feature_important)



######################### Predicting Hormone status of Breast Cancers #############################



######### ER (Estrogen receptor) (Positive or Negative) ##############

df_with_nine_dimensions <- data.frame(new_df_feature_important)

pca <- prcomp(df_with_nine_dimensions)

df_with_nine_dimensions <- pca$x[,1:9] #Selecting first 9 components
df_with_nine_dimensions <- data.frame(df_with_nine_dimensions)

# 1 means Positive ER  and 0 means Negative
df_with_nine_dimensions['ER'] <- ifelse(df$ER.Status == 'Positive', 1, 0)


#Creating a model by taking random samples multiple times
set.seed(123133456)
trainResult <- rep(0, 5)
testResult <- rep(0, 5)
for (i in c(1:5)) {
  sample_size <- sample.split(df_with_nine_dimensions, SplitRatio = 8 / 10)
  train <- subset(df_with_nine_dimensions, sample_size == T)
  test <- subset(df_with_nine_dimensions, sample_size == F)
  
  table(train$ER)
  
  logistic_model_train <-
    glm(formula =  ER ~ .,
        family = "binomial",
        data = train)
  summary(logistic_model_train)
  
  #Train set Predictions
  tainset_predicted_values <-
    predict(logistic_model_train, train, type = 'response')
  predicted <- ifelse(tainset_predicted_values > 0.5, 1, 0)
  confusion_matrix_train <-
    table(Predicted = predicted, Actual = train$ER)
  
  trainResult[i] <-
    sum(diag(confusion_matrix_train)) / sum(confusion_matrix_train)
  
  #Test Set Predictions
  testset_predicted_values <-
    predict(logistic_model_train, test, type = 'response')
  predicted <- ifelse(testset_predicted_values > 0.5, 1, 0)
  confusion_matrix_test <-
    table(Predicted = predicted, Actual = test$ER)
  
  testResult[i] <-
    sum(diag(confusion_matrix_test)) / sum(confusion_matrix_test)
  
}
Actual <- ifelse(test$ER == 0, 'ER Negative', 'ER Positive')
Predicted <- ifelse(predicted == 0, 'ER Negative', 'ER Positive')
table(Actual, Predicted)

cat("Average Train Accuracy = ",mean(trainResult))

cat("Average Test Accuracy = ",mean(testResult))



######### PR (Progestrone receptor) (Positive or Negative) ##############

df_with_nine_dimensions['ER'] <- NULL

df_with_nine_dimensions['PR'] <- ifelse(df$PR.Status == 'Positive', 1, 0)


df_with_nine_dimensions$PR%>% table()

trainResult <- rep(0, 5)
testResult <- rep(0, 5)
for (i in c(1:5)) {
  sample_size <- sample.split(df_with_nine_dimensions, SplitRatio = 9 / 10)
  train <- subset(df_with_nine_dimensions, sample_size == T)
  test <- subset(df_with_nine_dimensions, sample_size == F)
  
  
  logistic_model_train <-
    glm(formula =  PR ~ .,
        family = "binomial",
        data = train)
  summary(logistic_model_train)
  
  #Train set Predictions
  tainset_predicted_values <-
    predict(logistic_model_train, train, type = 'response')
  predicted <- ifelse(tainset_predicted_values > 0.5, 1, 0)
  confusion_matrix_train <-
    table(Predicted = predicted, Actual = train$PR)
  
  trainResult[i] <-
    sum(diag(confusion_matrix_train)) / sum(confusion_matrix_train)
  
  #Test Set Predictions
  testset_predicted_values <-
    predict(logistic_model_train, test, type = 'response')
  predicted <- ifelse(testset_predicted_values > 0.5, 1, 0)
  confusion_matrix_test <-
    table(Predicted = predicted, Actual = test$PR)
  
  testResult[i] <-
    sum(diag(confusion_matrix_test)) / sum(confusion_matrix_test)
  
}

Actual <- ifelse(test$PR == 0, 'PR Negative', 'PR Positive')
Predicted <- ifelse(predicted == 0, 'PR Negative', 'PR Positive')
cat("Average train Accuracy = ",mean(trainResult))
cat("Average Test Accuracy = ",mean(testResult))




################## HER2 #######################

df_with_nine_dimensions$PR <- NULL
df_with_nine_dimensions['HER2'] <- ifelse(df$HER2.Final.Status == 'Negative', 0, 1)

set.seed(123123)
trainResult <- rep(0, 5)
testResult <- rep(0, 5)

for (i in c(1:5)) {
  sample_size <- sample.split(df_with_nine_dimensions, SplitRatio = 8 / 10)
  train <- subset(df_with_nine_dimensions, sample_size == T)
  test <- subset(df_with_nine_dimensions, sample_size == F)
  
  logistic_model_train <-
    glm(formula =  HER2 ~ .,
        family = "binomial",
        data = train)
  summary(logistic_model_train)
  
  #Train set Predictions
  tainset_predicted_values <-
    predict(logistic_model_train, train, type = 'response')
  predicted <- ifelse(tainset_predicted_values > 0.5, 1, 0)
  confusion_matrix_train <-
    table(Predicted = predicted, Actual = train$HER2)
  
  trainResult[i] <-
    sum(diag(confusion_matrix_train)) / sum(confusion_matrix_train)
  
  #Test Set Predictions
  testset_predicted_values <-
    predict(logistic_model_train, test, type = 'response')
  predicted <- ifelse(testset_predicted_values > 0.5, 1, 0)
  confusion_matrix_test <-
    table(Predicted = predicted, Actual = test$HER2)
  
  testResult[i] <-
    sum(diag(confusion_matrix_test)) / sum(confusion_matrix_test)
  
}

Actual <- ifelse(test$HER2 == 0, 'HER2 Negative', 'HER2 Positive')
Predicted <- ifelse(predicted == 0, 'HER2 Negative', 'HER2 Positive')

cat("Average Train accuracy = ",mean(trainResult))


cat("Average Test accuracy = ",mean(testResult))




######## k-mean for ER, PR, HER2 positive #####

##### For K = 2 #########

set.seed(778513)
model <- kmeans(new_df_feature_important, centers = 2, iter.max = 500, nstart = 5)
model$betweenss/model$totss

fviz_cluster(model, new_df_feature_important)

table(model$cluster)

df_with_nine_dimensions %>% ggplot(aes(x = PC1, y = PC2)) + geom_point(color = "peru")


table(Clusters = model$cluster, df$ER) 
table(Clusters = model$cluster, df$PR.Status)
table(Clusters = model$cluster, df$HER2.Final.Status)



######## K-mean for Breast Cancer Categories for K = 4 ########

original_df <- new_df_feature_important
set.seed(778513)
model <- kmeans(original_df, centers = 4, iter.max = 500, nstart = 5)
model$betweenss/model$totss
fviz_cluster(model, original_df)
table(model$cluster, df$PAM50.mRNA)
treatement <- new_df_feature_important


####### Model to predict breast cancer group (SVM) ##########

library(e1071)

breast_cancer_patients <- new_df_feature_important
breast_cancer_patients <- data.frame(breast_cancer_patients)
breast_cancer_patients$cancerType <- df$PAM50.mRNA

#This data set contain data about proteins and breast cancer group
head(breast_cancer_patients)
table(breast_cancer_patients$cancerType)

#This model has overfiting problem
set.seed(778)
sample_size <-sample.split(breast_cancer_patients, SplitRatio = 8 / 10)
train <- subset(breast_cancer_patients, sample_size == T)
test <- subset(breast_cancer_patients, sample_size == F)

svmfit = svm(
  cancerType ~ .,
  data = train,
  kernel = "linear",
  cost = 1,
  scale = TRUE
)

#train predictions
predicted_values <- predict(svmfit, train)
confusion_matrix_test <-
  table(Predicted = predicted_values, Actual = train$cancerType)
print(confusion_matrix_test)
sum(diag(confusion_matrix_test)) / sum(confusion_matrix_test)

#test predictions
predicted_values <- predict(svmfit, test)
confusion_matrix_test <-
  table(Predicted = predicted_values, Actual = test$cancerType)
print(confusion_matrix_test)
sum(diag(confusion_matrix_test)) / sum(confusion_matrix_test)



##### SVM model for Comparison with k-mean  #######

breast_cancer_patients$cancerType <- as.character(breast_cancer_patients$cancerType)
breast_cancer_patients <- breast_cancer_patients[which(breast_cancer_patients$cancerType != 'HER2-enriched'),]
breast_cancer_patients$cancerType <- as.factor(breast_cancer_patients$cancerType)

pca <- prcomp(breast_cancer_patients[-length(breast_cancer_patients)])
features <- pca$x[,1:2]
features <- data.frame(features)
features$cancerType <- breast_cancer_patients$cancerType



svmfit <- svm(
  cancerType ~ .,
  data = features,
  kernel = "radial",
  cost = 1,
  scale = TRUE
)


predictions <- predict(svmfit, features)
confusion_matrix_test <-
  table(Predicted = predictions, Actual = features$cancerType)
print(confusion_matrix_test)
sum(diag(confusion_matrix_test)) / sum(confusion_matrix_test)



##### SVM model for making predictions ########
set.seed(99881)
sample_size <-sample.split(features, SplitRatio = 8 / 10)
train <- subset(features, sample_size == T)
test <- subset(features, sample_size == F)

svmfit <- svm(
  cancerType ~ .,
  data = train,
  kernel = "radial",
  cost = 1,
  scale = TRUE
)


plot(svmfit , train)

#train set predictions
predictions <- predict(svmfit, train)
confusion_matrix_test <-
  table(Predicted = predictions, Actual = train$cancerType)
print(confusion_matrix_test)
sum(diag(confusion_matrix_test)) / sum(confusion_matrix_test)

#test set predictions
predictions <- predict(svmfit, test)
confusion_matrix_test <-
  table(Predicted = predictions, Actual = test$cancerType)
print(confusion_matrix_test)

sum(diag(confusion_matrix_test)) / sum(confusion_matrix_test)


##################### K-mean (without HER2) #####################

fviz_nbclust(breast_cancer_patients[-length(breast_cancer_patients)], kmeans, method = "wss") +
  geom_vline(xintercept = 2, linetype = 2) + # add line for better visualization
  labs(subtitle = "Elbow method") # add subtitle

# Silhouette method
fviz_nbclust(breast_cancer_patients[-length(breast_cancer_patients)], kmeans, method = "silhouette") +
  labs(subtitle = "Silhouette method")


# NbClust function
nbclust_out <- NbClust(breast_cancer_patients[-length(breast_cancer_patients)],
                       distance = "euclidean",
                       min.nc = 2, max.nc = 4,
                       method = "kmeans")


model <- kmeans(breast_cancer_patients[-length(breast_cancer_patients)], centers = 3)
model$betweenss/model$totss
fviz_cluster(model, breast_cancer_patients[-length(breast_cancer_patients)])
table(model$cluster, breast_cancer_patients$cancerType) 







###################################################### Hierarchical ##########################################


######## Ward Method #######

#K = 4
dis_matrix <- dist(treatement, method = "euclidean")
d.ward <- hclust(dis_matrix, method = "ward.D2")
plot(d.ward, main = "Ward Linkage Clustering")
res <- cutree(d.ward, k = 4)
table(Clusters = res, df$PAM50.mRNA)

#K = 3 and without HER2
dis_matrix <- dist(breast_cancer_patients[-length(breast_cancer_patients)], method = "euclidean")
d.ward <- hclust(dis_matrix, method = "ward.D2")
plot(d.ward, main = "Ward Linkage Clustering")
res <- cutree(d.ward, k = 3)
table(Clusters = res, breast_cancer_patients$cancerType)



##### Single Linkage #####
dis_matrix <- dist(treatement, method = "euclidean")
d.single <- hclust(dis_matrix, method = "single")
plot(d.single, main = "Single Linkage Clustering")
res <- cutree(d.single, k = 4)
table(Clusters = res, df$PAM50.mRNA)


##### Average Linkage #####
dis_matrix <- dist(treatement, method = "euclidean")
d.average <- hclust(dis_matrix, method = "average")
plot(d.average, main = "Average Linkage Clustering")
res <- cutree(d.average, k = 4)
table(Clusters = res, df$PAM50.mRNA)


###### Complete Linkage ###########
dis_matrix <- dist(treatement, method = "euclidean")
d.complete <- hclust(dis_matrix, method = "complete")
plot(d.complete, main = "Complete Linkage Clustering")
res <- cutree(d.complete, k = 4)
table(Clusters = res, df$PAM50.mRNA)




