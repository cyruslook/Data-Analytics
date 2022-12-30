library(ggplot2)
library(dplyr)
library(corrplot)
library(lmtest)
library('aod')
library(generalhoslem)
library(caTools)

loan_data <- read.csv('./Documents/PredictAnalytics/GroupProjects/loan_data_set.csv')

head(loan_data)
loan_data$Loan_ID <- NULL
#Number of missing values in each column
NAcol <- which(colSums(is.na(loan_data)) > 0);NAcol
sort(colSums(sapply(loan_data[NAcol], is.na)), decreasing = TRUE)


#Getting numerical variables to check the coorelation
numericVars <- which(sapply(loan_data, is.numeric))
numericVarNames <- names(numericVars)

all_numVar <- loan_data[, c(numericVarNames)]
M <- cor(all_numVar, use="pairwise.complete.obs")

#Corrplot suggested no auto correlation
corrplot(M, method = "number", type = "upper")

#Total 50 NULL values in the credit history
loan_data %>% ggplot(aes(x = as.factor(Credit_History), fill = as.factor(Credit_History))) + geom_bar(stat = 'count') + geom_label(stat='count', aes(label=..count..))

#We can see from the above plot that most of the peopel in the dataset have credit history, now less check how many of the has approved loan status
loan_data %>% ggplot(aes(x = as.factor(Credit_History), fill = Loan_Status)) + 
  geom_bar(stat='count', position='dodge') + theme_grey() + geom_label(stat='count', aes(label=..count..))

#As we can see for most of the people who have credit history, their loan is approved, 
#whereas among the people with no credit history, only 7 out of 89 have thier loan approved

#Going to replace the NA values in credit_history accordin to the loan_status. i.e if loan status is 1 means credit_history = 1
table(loan_data$Credit_History)
loan_data$Credit_History <- ifelse(loan_data$Credit_History %in% c(0, 1), 
                                   loan_data$Credit_History,
                                   ifelse(loan_data$Loan_Status == 'Y', 1, 0))



#Check loan approval according to their education level
loan_data %>% ggplot(aes(x = as.factor(Education), fill = Loan_Status)) + 
  geom_bar(stat='count', position='dodge') + theme_grey() + geom_label(stat='count', aes(label=..count..))

#Check loan approval according to their Enployment Status
table(loan_data$Self_Employed)
loan_data %>% ggplot(aes(x = as.factor(Self_Employed), fill = Loan_Status)) + 
  geom_bar(stat='count', position='dodge') + theme_grey() + geom_label(stat='count', aes(label=..count..))


#Loan Amount Term
table(loan_data$Loan_Amount_Term)

#LoanAmount Term Vs LoanAmount

#Also we can see that there is no observation in the dataset for which both the loan_amount_term and loanAmount is NULL
#Means if there is clear relationship between loan amount and loan_amount_term, we can easily fill the missing values


#Property area vs LoanAmount
loan_data %>% ggplot(aes(x = Property_Area, y = LoanAmount, color = Property_Area)) + geom_boxplot() + geom_jitter()
table(loan_data$Loan_Amount_Term)



################### T - Test #####################

########Loan Amount##############
res.ftest <- var.test(LoanAmount ~ Loan_Status, data = loan_data[!is.na(loan_data$LoanAmount),])
res.ftest
loan_data %>% ggplot(aes(x = Loan_Status, y = LoanAmount, color = Loan_Status)) + geom_boxplot() + geom_jitter()

res <- t.test(LoanAmount ~ Loan_Status, data = loan_data[!is.na(loan_data$LoanAmount),], var.equal = TRUE)
res
###### Applicant Income #########
res.ftest <- var.test(ApplicantIncome ~ Loan_Status, data = loan_data[!is.na(loan_data$LoanAmount),])
res.ftest

loan_data %>% ggplot(aes(x = Loan_Status, y = ApplicantIncome, color = Loan_Status)) + geom_boxplot() + geom_jitter()

res <- t.test(ApplicantIncome ~ Loan_Status, data = loan_data[!is.na(loan_data$LoanAmount),], var.equal = FALSE)
res
###### CoapplicantIncome ############
res.ftest <- var.test(CoapplicantIncome ~ Loan_Status, data = loan_data[!is.na(loan_data$LoanAmount),])
res.ftest

loan_data %>% ggplot(aes(x = Loan_Status, y = CoapplicantIncome, color = Loan_Status)) + geom_boxplot() 

res <- t.test(CoapplicantIncome ~ Loan_Status, data = loan_data[!is.na(loan_data$LoanAmount),], var.equal = FALSE)
res




loan_data[!is.na(loan_data$LoanAmount),] %>%
  ggplot(aes(x = Loan_Status, y = LoanAmount, color = Loan_Status)) + geom_boxplot() + geom_jitter(alpha = 0.8)
#Although there is no clear relation of LoanAmount with any other variable, we are going to replace missing values in LoanAmount according to their median of loan status, so that we will not be having any issue of outlier


loan_data <- loan_data %>%
  group_by(Loan_Status) %>%
  mutate(LoanAmount = ifelse(
    is.na(LoanAmount) ,
    median(LoanAmount, na.rm = TRUE),
    LoanAmount
  ))


#Loan_amount_term according to loan_status
loan_data %>% ggplot(aes(x = as.factor(Loan_Amount_Term), fill = Loan_Status)) + 
  geom_bar(stat='count', position='dodge') + theme_grey() + geom_label(stat='count', aes(label=..count..))


#just like loan amount, there is no clear relationship of loan_amount_term with other variables
#As we have considered the Loan_Amount_Term as a discrete variable, we will replace missing values in the Loan_Amount_Term with maximum repeated value for the same column
get_maximum_repeted_value <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

loan_data$Loan_Amount_Term <-ifelse(is.na(loan_data$Loan_Amount_Term), get_maximum_repeted_value(loan_data$Loan_Amount_Term), loan_data$Loan_Amount_Term) 
table(loan_data$Loan_Amount_Term)

#Above table suggested that the data for the loan_amount_term values(12, 36, 60, 84, 120, 240, 300, 480) is not enough to train our model.
# These rows like like a outlier in our dataset. So will remove these rows
loan_data <- loan_data[loan_data$Loan_Amount_Term == 180 | loan_data$Loan_Amount_Term == 360,]


options(dplyr.width = Inf) ##Setting width infinite, so that we can all the columns
head(loan_data)

table(loan_data$Gender)
table(loan_data$Married)
table(loan_data$Dependents)
table(loan_data$Education)
table(loan_data$Loan_Amount_Term)
table(loan_data$Property_Area)
table(loan_data$Self_Employed)

#There are 3 missing values in Married status and 12 in Dependents. Moreover, missing values in
#Married status belong to 1 of the 12 missing observations in Dependents
#Na in Dependents means no dependent & missing in married status means Not married

loan_data %>% ggplot(aes(x = Dependents, fill = Dependents)) + geom_bar() + theme_grey() + geom_label(stat='count', aes(label=..count..))

loan_data$Dependents <- as.character(loan_data$Dependents)
loan_data$Dependents <- ifelse(loan_data$Dependents == "", "0", loan_data$Dependents)


loan_data$Married <- as.character(loan_data$Married)
loan_data$Married <- ifelse(loan_data$Married == "", "No", loan_data$Married)


##Replacing the missing values for gender column
loan_data$Gender <- as.character(loan_data$Gender)
loan_data[loan_data$Gender != '',] %>%
  ggplot(aes(x = Gender, fill = Loan_Status)) + geom_bar(stat='count', position='dodge') + theme_grey() + geom_label(stat='count', aes(label=..count..))

#It is clear from the plot that loan approval is not dependent on the gender(trend is same for both male and Female)
#As most of the Gender are Male, we are going to replace missing values in gender with male

loan_data$Gender <- ifelse(loan_data$Gender == "", "Male", loan_data$Gender)
chisq.test(loan_data$Gender, loan_data$Loan_Status)


#Replacing the missing values in self employed
loan_data$Self_Employed <- as.character(loan_data$Self_Employed)

loan_data[loan_data$Self_Employed != '',] %>%
  ggplot(aes(x = Self_Employed, fill = Loan_Status)) + geom_bar(stat='count', position='dodge') + theme_grey() + geom_label(stat='count', aes(label=..count..))

#It is clear from the plot that Self_Employed is not dependent on the loan Approval, below chi-square is proof to this.As most of the customers are Not self employed, we will replace the missing values with "No"
loan_data$Loan_Status <- as.character(loan_data$Loan_Status)
chisq.test(loan_data$Self_Employed, loan_data$Loan_Status)

loan_data$Self_Employed <- ifelse(loan_data$Self_Employed == "", "No", loan_data$Self_Employed)


#By this time all the missing values are filled
sum(is.na(loan_data))

dim(loan_data)

#Converting some character variables back to factor
loan_data$Loan_Status <- as.factor(loan_data$Loan_Status)
loan_data$Self_Employed <- as.factor(loan_data$Self_Employed)
loan_data$Married <- as.factor(loan_data$Married)
loan_data$Dependents <- as.factor(loan_data$Dependents)
loan_data$Gender <- as.factor(loan_data$Gender)


loan_data <- data.frame(loan_data)
set.seed(101)
sample_size <- sample.split(loan_data, SplitRatio = 8/10)
train <-subset(loan_data, sample_size == T)
test <-subset(loan_data, sample_size == F)


#Variable Selection

Intercept_only_model <-glm(Loan_Status~1,data=train,family = "binomial")

full_model <- train %>% glm(formula =  Loan_Status ~ ., family = "binomial")


step_wise_selection <- step(Intercept_only_model, direction = 'both', scope = formula(full_model))

summary(step_wise_selection)

#Only Credit_History & Property_Area are significant predictors according to step_wise
final_model <- glm(Loan_Status~Credit_History + Property_Area, data=train, family = "binomial")
summary(final_model)

########Interaction Plot#########

interaction.plot(x.factor = train$Credit_History,    # variable to plot on x-axis
                 trace.factor = train$Property_Area, # variable to specify "traces"; here, lines
                 response = as.integer(train$Loan_Status),    # variable to plot on y-axis
                 fun = mean,  # summary statistic to be plotted for response variable
                 type = "l",     # type of plot, here "l" for lines
                 ylab = "",
                 xlab = "",
                 col = c("blue4", "red4"),
                 lty = 1,  # line type
                 lwd = 2,  # line width
                 trace.label = "",  # label for legend
                 xpd = FALSE) #,  # 'clip' legend at border




#There is an interaction between gender and education

final_model_interaction_terms <- train %>% glm(formula =  Loan_Status ~ Credit_History + Property_Area + Credit_History*Property_Area, family = "binomial")
summary(final_model_interaction_terms)

#-------Wald Test-------------------
wald.test(b=coef(final_model),Sigma = vcov(final_model),Terms = 2)
wald.test(b=coef(final_model),Sigma = vcov(final_model),Terms = 3)
wald.test(b=coef(final_model),Sigma = vcov(final_model),Terms = 4)


#Likelihood Ratio Test for two different model
pchisq(388.36 - 387.08, 424 - 422, lower.tail = F)



#ROC
library(pROC)
roc_curve <- roc(train$Loan_Status, 
    final_model$fitted.values, 
    plot = TRUE, 
    legacy.axes = T, 
    percent = T, 
    lwd = 5, 
    col = 'blue',
    xlab = "False Positive Percentage",
    ylab = "True Positive Percentage",
    print.auc = T)

plot.roc(train$Loan_Status, 
         final_model_interaction_terms$fitted.values, 
         percent = T,
         col = 'green',
         print.auc.y = 40,
         add = T,
         print.auc = T,
         )


roc_data <- data.frame(
  true_positive = roc_curve$sensitivities,
  false_positive = (100 - roc_curve$specificities),
  thresholds = roc_curve$thresholds
)

roc_data

#Model Performance
tainset_predicted_values <- predict(final_model, train, type='response')
y_pred_train = ifelse(tainset_predicted_values > 0.77, 1, 0)

confusion_matrix_train <- table(Predicted = y_pred_train, Actual=train$Loan_Status)
sum(diag(confusion_matrix_train))/sum(confusion_matrix_train)

sensitivity <- confusion_matrix_train[2,2]/(confusion_matrix_train[2,2] + confusion_matrix_train[1,2]);sensitivity

specificity <- confusion_matrix_train[1,1]/(confusion_matrix_train[1,1] + confusion_matrix_train[2,1]);specificity


testset_predicted_values <- predict(final_model, test, type='response')
y_pred_test = ifelse(testset_predicted_values > 0.77, 1, 0)

confusion_matrix_test <- table(Predicted = y_pred_test, Actual=test$Loan_Status)

sensitivity <- confusion_matrix_test[2,2]/(confusion_matrix_test[2,2] + confusion_matrix_test[1,2]);sensitivity

specificity <- confusion_matrix_test[1,1]/(confusion_matrix_test[1,1] + confusion_matrix_test[2,1]);specificity

#Sensitivity = TP/(TP + FN)
#Specificity = TN/(TN + FP)

Accuracy <- sum(diag(confusion_matrix_test))/sum(confusion_matrix_test);Accuracy



#Confidence interval for parameters
confint(final_model)




#Hosmer and Lemeshow test
final_model_train <- glm(Loan_Status~Credit_History + Property_Area, data=train, family = "binomial")
final_model_test <- glm(Loan_Status~Credit_History + Property_Area, data=test, family = "binomial")
summary(final_model_train)
summary(final_model_test)

logitgof(train$Loan_Status, fitted(final_model_train), g=10, ord=FALSE)
logitgof(test$Loan_Status, fitted(final_model_test), g=10, ord=FALSE)

#https://stats.stackexchange.com/questions/169438/evaluating-logistic-regression-and-interpretation-of-hosmer-lemeshow-goodness-of
#https://www.youtube.com/watch?v=_s3XNX3TlZM



