# primary functions
library(tibble)
library(dplyr)
library(FactoMineR)
library(caret)
library(plm)

#' @param df dataframe (decide if rows have been randomized)
#'    randomized outside of function if mix R/python methods
#' @param K number of folds for cross-validation
#' @param nperfmeas number of performance measures used to compare methods
#' @param seed random number of seed for reproducibility
#' @return matrix of performance measures by fold, vector of averages
crossValidationCont = function(df,K,nperfmeas=6,seed)
{ set.seed(seed)
  # apply transformations
  df <- transform(df)
  
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
  commdf <- ClassSplit(df)[[1]]%>% medianImp()
  resdf <- ClassSplit(df)[[2]]%>% medianImp()
  
  # commercial
  commtrain = commdf[-indices,]
  commholdout = commdf[indices,]
  
  #PredfestcomModel <- feols(site_eui~minavg+maxavg+energy_star_rating+floor_area+facility_type+Year_Factor+State_Factor+cooling_degree_days+heating_degree_days+precipitation_inches+days_above_110F+days_below_30F,commtrain)
  PredfestcomModel <- causalmodel(commtrain)
  
  restrain = resdf[-indices,] 
  resholdout = resdf[indices,] 
  #PredfestresModel <- feols(site_eui~minavg+maxavg+energy_star_rating+facility_type+Year_Factor+State_Factor+cooling_degree_days+heating_degree_days+precipitation_inches+days_above_110F+days_below_30F,restrain)
  PredfestresModel <- causalmodel(restrain)
  #pred = myPredict(obj, newdata=holdout)  # could include 50% and 80% intervals
  
  comperf50meas[k,] = probCheck(PredfestcomModel,commholdout,0.50)$summary
  comperf80meas[k,] = probCheck(PredfestcomModel,commholdout,0.80)$summary
  resperf50meas[k,] = probCheck(PredfestresModel,resholdout,0.50)$summary
  resperf80meas[k,] = probCheck(PredfestresModel,resholdout,0.80)$summary
  #perfMeas(pred, holdout$y)
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

# functions

#' Causual model implementation
#'
#' @param df 
#'
#' @return casual model object
#' @export
#'
#' @examples
causalmodel <- function(df){
  X <- df %>% select(-c('id'))#,'Year_Factor','State_Factor','facility_type'))#%>% mutate_if(is.character,as.numeric) %>% mutate_all(~. -mean(.,na.rm=TRUE))
  #X$State_Factor = df$State_Factor
  #X$Year_Factor = df$Year_Factor
  #X$facility_type = df$facility_type
  
  # auxillary
  Xaux <- lm(avg_temp~cooling_degree_days*days_above_110F+heating_degree_days*days_below_30F+precipitation_inches,X)
  X$avg_temp <- Xaux$fitted.values
  X <- X %>% select(-c("cooling_degree_days","days_above_110F","heating_degree_days","days_below_30F","precipitation_inches"))
  
  Xmod <- lm(site_eui~minavg+maxavg+energy_star_rating+facility_type+floor_area+ELEVATION+State_Factor+Year_Factor+avg_temp+floor_area*ELEVATION+floor_area*energy_star_rating,X)
  return(Xmod)
}


#' Binning State Factor
#'
#' @param df 
#' @param ntrain 
#'
#' @return transformed State factors in 3 separate bins H, M, C. According to climates.
#' @export
#'
#' @examples
statebins <- function(df,ntrain){
  weather= rep("M",ntrain)
  weather[df$State_Factor == 'State_4']='C';weather[df$State_Factor == 'State_6']='C';weather[df$State_Factor == 'State_1']='H';weather[df$State_Factor == 'State_11']='H';
  weather = as.factor(weather); df$State_Factor = weather
  return(df)
}

#' No outlier transformation using IQR
#'
#' @param df 
#'
#' @return subset of the dataset after removing outliers
#' @export
#'
#' @examples
noout <- function(df){
  iqr = IQR(df$site_eui)
  Q1 <- quantile(df$site_eui,.25)
  Q3 <- quantile(df$site_eui,.75)
  res <- subset(df, df$site_eui > (Q1 - 1.5*iqr) & df < (Q3 + 1.5*iqr))
  return(na.omit(res))
}

#' TRANSFORMATION 
#'
#' @param df 
#' @param subset 
#'
#' @return df with all the transformations in place.
#' @export
#'
#' @examples
transform <- function(df,subset = TRUE){
  if (subset==TRUE){
    df <- df
    ntrain <- nrow(df)/2
    df <- df[sample(nrow(df),ntrain),]
    df <-  statebins(df,ntrain)
  }
  # transformations
  
  df$floor_area <- log(df$floor_area)
  # remove outliers
  #df$site_eui <- log(df$site_eui)
  # using IQR
  df <- noout(df)
  # days_data
  daysDatNames <- df %>% select(contains('days')) %>% names()
  
  # creating auxillarydf
  #auxDf <- df %>% select(c(daysDatNames,'avg_temp',"cooling_degree_days","heating_degree_days","precipitation_inches","snowfall_inches","snowdepth_inches","building_class"))
  
  # removing missing values
  
  #df <- df %>% select(-c('id','direction_peak_wind_speed' ,'max_wind_speed' ,'days_with_fog',"cooling_degree_days","heating_degree_days")) 
  # removed percip,snow,days
  
  #df <- df %>% select(-names(df)[45:57])
  # min and max names
  minNames <- df %>% select(contains("min")) %>% names()
  maxNames <- df%>% select(contains("max")) %>% names()
  avgNames <- df%>% select(contains("_avg_")) %>% names()
  
  # create a minimum and maximum averagese.
  df <- df %>% mutate(minavg = (january_min_temp+february_min_temp+march_min_temp+october_min_temp+november_min_temp+december_min_temp)/4)
  df <- df %>% mutate(maxavg = (april_max_temp+may_max_temp+june_max_temp+july_max_temp+august_max_temp+september_max_temp)/4)
  
  #return(list(df,auxDf))
  return(df)
  
}

#' Median Imputations
#'
#' @param df 
#'
#' @return df with year built and energy star rating na values converted to median values
#' @export
#'
#' @examples
medianImp <- function(df){
  df$year_built[is.na(df$year_built)]= summary(df$year_built)[[3]]
  df$energy_star_rating[is.na(df$energy_star_rating)] = summary(df$energy_star_rating)[[3]]
  return(df)
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
  commdf <- df %>% filter(building_class=='Commercial')%>% select(-'building_class')
  resDf <- df %>% filter(building_class=="Residential") %>% select(-'building_class')
  return(list(commdf,resDf))
}

#' Density plots
#'
#' @param df 
#' @param model 
#'
#' @return plots site_eui(red) and estimated site_eui density, to compare model performance
#' @export
#'
#' @examples
plot_results <- function(df,model){
  preds <- predict(model,df)
  plot(density(preds));
  lines(density(df$site_eui),col=2,type='h',main='predictions vs actual');
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