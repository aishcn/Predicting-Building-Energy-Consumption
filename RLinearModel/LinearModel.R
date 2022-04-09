
# Fri Apr  8 21:22:28 PDT 2022

# loading the libraries
source('Functions1.R')

# import the dataset
bigdf <- read.csv('2_df_transform.csv')
print(summary(causalmodel2(bigdf)))

# creating data partitions
dfintrain <- createDataPartition(bigdf$site_eui,p=.8,list=FALSE)

trainDf <- bigdf[dfintrain,]
testDf <- bigdf[-dfintrain,]

# fit control
fitControl <- trainControl(## 10-fold CV
                           method = "repeatedcv",
                           number = 5,
                           ## repeated ten times
                           repeats = 3,
)

model <- train(site_eui~SummerTemp+WinterTemp+FallTemp+energy_star_rating+facility_type+floor_area+as.factor(Year_Factor)+avg_temp+floor_areaxELEVATION+floor_areaxenergy_star_rating, data = trainDf, 
                 method = "lm", 
                 preProcess = c('center', 'scale'),
                  trControl = fitControl,
                          metric = c('RMSE,MA,Rsquared'))
getTrainPerf(model) %>% print()

# cross validation that outputs interval score.
crossvalResults = crossValidationCont(bigdf,5,nperfmeas = 4,seed=12)

print(crossvalResults)

fiftyAvgPred = (crossvalResults$resavgperf50meas*0.57) + (crossvalResults$comavgperf50meas*0.43)
print('50 % prediction interval for model comparison')
print(fiftyAvgPred)
eightyAvgPred = (crossvalResults$resavgperf80meas*0.57) + (crossvalResults$comavgperf80meas*0.43)
print('80 % prediction interval for model comparison')
print(eightyAvgPred)
