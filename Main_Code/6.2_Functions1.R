# primary functions
library(dplyr)
library(caret)

#' @param df dataframe (decide if rows have been randomized)
#'    randomized outside of function if mix R/python methods
#' @param K number of folds for cross-validation
#' @param nperfmeas number of performance measures used to compare methods
#' @param seed random number of seed for reproducibility
#' @return matrix of performance measures by fold, vector of averages
crossValidationCont = function(df,K,nperfmeas=6,seed)
{ set.seed(seed)
  
  n = nrow(df)
  nhold = round(n/K) # size of holdout set 
  iperm = sample(n)
  comperf50meas = matrix(0,K,nperfmeas)
  comperf80meas = matrix(0,K,nperfmeas)
  resperf50meas = matrix(0,K,nperfmeas)
  resperf80meas = matrix(0,K,nperfmeas)
  
  for(k in 1:K)
  { indices = (((k-1)*nhold+1):(k*nhold))
  if( k==K ) indices = (((k-1)*nhold+1):n)
  indices = iperm[indices]
  
  # split between commercial and residential 
  commdf <- ClassSplit(df)[[1]]
  resdf <- ClassSplit(df)[[2]]
  
  # commercial
  commtrain = commdf[-indices,]
  commholdout = commdf[indices,]
  
  PredfestcomModel <- causalmodel2(commtrain)
  
  # residential
  restrain = resdf[-indices,] 
  resholdout = resdf[indices,] 
  PredfestresModel <- causalmodel2(restrain)
  
  # Commercial probability measures
  comperf50meas[k,] = probCheck(PredfestcomModel,commholdout,0.50)$summary
  comperf80meas[k,] = probCheck(PredfestcomModel,commholdout,0.80)$summary
  # Residential probability measures
  resperf50meas[k,] = probCheck(PredfestresModel,resholdout,0.50)$summary
  resperf80meas[k,] = probCheck(PredfestresModel,resholdout,0.80)$summary

  }
  # 50
  comavgperf50meas = apply(comperf50meas,2,mean)
  resavgperf50meas = apply(resperf50meas,2,mean)
  # 80 
  comavgperf80meas = apply(comperf80meas,2,mean)
  resavgperf80meas = apply(resperf80meas,2,mean)
  
  list(comperfmeas50byfold=comperf50meas,resperfmeas50byfold = resperf50meas, comavgperf50meas=comavgperf50meas,resavgperf50meas = resavgperf50meas,comperfmeas80byfold=comperf80meas,resperfmeas80byfold = resperf80meas, comavgperf80meas=comavgperf80meas,resavgperf80meas = resavgperf80meas)
}


# prob matrix
#' Computes probability matrix check
#'
#' @param model 
#' @param df 
#' @param level 
#'
#' @return interval score
#' @export
#'
#' @examples
probCheck <- function(model,df,level){
  problist <- predict(model,df,interval='prediction',level=level) %>% na.omit()
  return(intervalScore(problist,na.omit(df$site_eui),level))
}


#' Causal Model: as per equations presented in the report.
#'
#' @param df 
#'
#' @return model object
#' @export
#'
#' @examples
causalmodel2 <- function(df){
  
  # auxillary created an estimator for the average temperature.
  Xaux <- lm(avg_temp~freezing_days+warm_days*cold_days+snow_rain_inches,df)
  df$avg_temp <- Xaux$fitted.values
  
  Xmod <- lm(site_eui~SummerTemp+WinterTemp+FallTemp+energy_star_rating+facility_type+floor_area+as.factor(Year_Factor)+avg_temp+floor_areaxELEVATION+floor_areaxenergy_star_rating,df)
  return(Xmod)
}

#' Class Split
#'
#' @param df 
#'
#' @return two dataframes as a list commercial dataframe and residential dataframe. 
#' @export
#'
#' @examples
ClassSplit <- function(df){
  #commdf <- df %>% filter(building_class=='Commercial')%>% select(-'building_class')
  commdf <- df %>% filter(building_class==0)%>% select(-'building_class')
  #resDf <- df %>% filter(building_class=="Residential") %>% select(-'building_class')
  resDf <- df %>% filter(building_class==1) %>% select(-'building_class')
  return(list(commdf,resDf))
}


#' Interval Score
#'
#' @param predObj 
#' @param actual 
#' @param level 
#'
#' @return Interval Score as written by Harry Joe.
#' @export
#'
#' @examples
intervalScore = function(predObj,actual,level)
{ n = nrow(predObj)
alpha = 1-level
ilow = (actual<predObj[,2])  # underestimation
ihigh = (actual>predObj[,3]) # overestimation
sumlength = sum(predObj[,3]-abs(predObj[,2])) # sum of lengths of prediction intervals
sumlow = sum(predObj[ilow,2]-actual[ilow])*2/alpha
sumhigh = sum(actual[ihigh]-predObj[ihigh,3])*2/alpha
avglength = sumlength/n
IS = (sumlength+sumlow+sumhigh)/n # average length + average under/over penalties
cover = mean(actual>= predObj[,2] & actual<=predObj[,3])
summ = c(level,avglength,IS,cover) 
# summary with level, average length, interval score, coverage rate
imiss = which(ilow | ihigh)
list(summary=summ, imiss=imiss)
}
