# STAT 447B Group Project
## Predicting Building Energy Consumption
**Contributors**: ***Anjali Chauhan, Idris Hedayat, Sameer Shankar, Sumit Meghlani***

## Introduction
Climate change is a globally relevant, urgent, and multi-faceted issue heavily impacted by energy policy and infrastructure. Addressing climate change involves mitigation (i.e. mitigating greenhouse gas emissions) and adaptation (i.e. preparing for unavoidable consequences). Mitigation of GHG emissions requires changes to electricity systems, transportation, buildings, industry, and land use.

According to a report issued by the International Energy Agency (IEA), the lifecycle of buildings from construction to demolition were responsible for 37% of global energy-related and process-related CO2 emissions in 2020. Yet it is possible to drastically reduce the energy consumption of buildings by a combination of easy-to-implement fixes and state-of-the-art strategies. For example, retrofitted buildings can reduce heating and cooling energy requirements by 50-90 percent. Many of these energy efficiency measures also result in overall cost savings and yield other benefits, such as cleaner air for occupants. This potential can be achieved while maintaining the services that buildings provide.

**Dataset**: 

Created in collaboration with Climate Change AI (CCAI) and Lawrence Berkeley National Laboratory (Berkeley Lab).It consists of variables that describe building characteristics and climate and weather variables for the regions in which the buildings are located. Accurate predictions of energy consumption can help policymakers target retrofitting efforts to maximize emissions reductions.

- Train data_set has 75757 rows and 64 columns
- Test dataset has 9705 rows and 63 columns
- Only `State_Factor`, `building_class` and `facility_type` are categorical variables. Rest of the variables are all numerical.

**Description of variables**:

`id`: building id

`Year_Factor`: anonymized year in which the weather and energy usage factors were observed

`State_Factor`: anonymized state in which the building is located

`building_class`: building classification

`facility_type`: building usage type

`floor_area`: floor area (in square feet) of the building

`year_built`: year in which the building was constructed

`energy_star_rating`: the energy star rating of the building

`ELEVATION`: elevation of the building location

`january_min_temp`: minimum temperature in January (in Fahrenheit) at the location of the building

`january_avg_temp`: average temperature in January (in Fahrenheit) at the location of the building

`january_max_temp`: maximum temperature in January (in Fahrenheit) at the location of the building

`cooling_degree_days`: cooling degree day for a given day is the number of degrees where the daily average temperature exceeds 65 degrees Fahrenheit. Each month is summed to produce an annual total at the location of the building

`heating_degree_days`: heating degree day for a given day is the number of degrees where the daily average temperature falls under 65 degrees Fahrenheit. Each month is summed to produce an annual total at the location of the building

`precipitation_inches`: annual precipitation in inches at the location of the building

`snowfall_inches`: annual snowfall in inches at the location of the building

`snowdepth_inches`: annual snow depth in inches at the location of the building

`avg_temp`: average temperature over a year at the location of the building

`days_below_30F`: total number of days below 30 degrees Fahrenheit at the location of the building

`days_below_20F`: total number of days below 20 degrees Fahrenheit at the location of the building

`days_below_10F`: total number of days below 10 degrees Fahrenheit at the location of the building

`days_below_0F`: total number of days below 0 degrees Fahrenheit at the location of the building

`days_above_80F`: total number of days above 80 degrees Fahrenheit at the location of the building

`days_above_90F`: total number of days above 90 degrees Fahrenheit at the location of the building

`days_above_100F`: total number of days above 100 degrees Fahrenheit at the location of the building

`days_above_110F`: total number of days above 110 degrees Fahrenheit at the location of the building

`direction_max_wind_speed`: wind direction for maximum wind speed at the location of the building. Given in 360-degree compass point directions (e.g. 360 = north, 180 = south, etc.)

`direction_peak_wind_speed`: wind direction for peak wind gust speed at the location of the building. Given in 360-degree compass point directions (e.g. 360 = north, 180 = south, etc.)

`max_wind_speed`: maximum wind speed at the location of the building

`days_with_fog`: number of days with fog at the location of the building

**Target variable**:

`site_eui`: Target Site Energy Usage Intensity is the amount of heat and electricity consumed by a building as reflected in utility bills

**Task**: 

Analyze differences in building energy efficiency, creating models to predict building energy consumption. 

**Evaluation Metric**: 

The evaluation metric for this competition is Root Mean Squared Error (RMSE). The RMSE is commonly used measure of the differences between predicted values provided by a model and the actual observed values.

## Exploratory Data Analysis (EDA)

## Feature Engineering

## Model Training

## Model Evaluation

## Hyperparameter Tuning

## Prediction

## Conclusion
