###########################################################-
# Objective: Quantify influence of catchment characteristics
# on flood generating processes
# Author: Lina Stein, University of Bristol
# note: 
###########################################################-


# Packages ----------------------------------------------------------------
library(quantmod)
library(chron)
library(data.table)
library(plyr)
library(dplyr)

#Visualisation packages
library(ggplot2)
library(maps)
library(reshape)
library(scales)
library(ggpubr)
library(stringr)
library(gtable)
library(lemon)
library(ggcorrplot)
library(patchwork)

#Random forest and accumulated local effects
library(randomForest)
library(hydroGOF)
library(gtools)
library(parallel)
library(iml)
library(devtools)


plot_path =  'C:/Users/ls16959/OneDrive - University of Bristol/Documents/8_plots/CAMELS_processes/'

#Location of CAMELS files
CAMELS_path = 'C:/Users/ls16959/Data/streamflow_countries/US/Camels/'
#Location of additional files (elongation ratio, AWS)
path_additional = "C:/Users/ls16959/Data/02_CatchmentCharacteristics/"
#Save location for large output files
save_path = "C:/Users/ls16959/Data/04_Calculations/CAMELS/"

# Functions ---------------------------------------------------------------

#read flood classification functions (from Stein et al, 2019)
source("C:/Users/ls16959/OneDrive - University of Bristol/Documents/1_Data/Code_by_me/variable_method_flood_mech_functions.R")

#read in observed streamflow from CAMELS data
read_qobs = function(gauge_id, filepaths){
  #read observed Q and remove irregular empty space
  qobs1 <- readLines(qfiles[grep(gauge_id, filepaths)][1]) #same file in different folders, only use 1st
  qobs2 = strsplit(qobs1, "\t| ")
  qobs3 = lapply(qobs2, FUN = function(x){
    x[which(x!="")]
  })
  qobs4 = do.call(rbind, qobs3)
  qobs5 = qobs4[-1,] #remove head
  qobs_date = as.Date(apply(qobs5[,1:3],1, FUN = function(x){paste(x[1], x[2], x[3], sep = "-")}), format = "%Y-%m-%d") #create date column
  qobs_Q = as.numeric(qobs5[,12]) #get observed Q values
  qobs = data.frame(qobs_date, qobs_Q)
}

#read in climate data from CAMELS data
read_climate = function(gauge_id, filepaths){
  #read observed Q and remove irregular empty space
  qobs1 <- readLines(qfiles[grep(gauge_id, filepaths)][1]) #same file in different folders, only use 1st
  qobs2 = strsplit(qobs1, "\t| ")
  qobs3 = lapply(qobs2, FUN = function(x){
    x[which(x!="")]
  })
  qobs4 = do.call(rbind, qobs3)
  qobs5 = qobs4[-1,] #remove head
  colnames(qobs5) = qobs4[1,]
  
  P_CAMELS = data.frame(as.numeric(as.character(qobs5[,"PRCP"])))
  ET_CAMELS = data.frame(as.numeric(as.character(qobs5[,"ET"])))
  T_CAMELS = data.frame(as.numeric(as.character(qobs5[,"TAIR"])))
  
  climatedf = data.frame(P_CAMELS, ET_CAMELS, T_CAMELS)
  colnames(climatedf) = c('RAIM', 'ET', 'TAIR')
  return(climatedf)
}


#Find peaks over threshold
POT_func = function(qobsdf, POT_thresh, eventsperyear){
  #1. set low threshold (median?)
  #2. Find Peaks
  #check independence
  #3. determine average time to rise for five clean (not multi-peak events (for example seperate years or AMAX?))
  ###Input
  #qobsdf: discharge timeseries with date and flow
  #POT_thresh: threshold for POT analysis
  #eventsperyear: limit for in average how many events occur per year
  
  
  Qobsdf = qobsdf
  colnames(Qobsdf) = c("Qobs_date", "Qobs")
  peak_vec = findPeaks(Qobsdf[,2])-1 #is always off by one
  #find peaks with a large difference to next peak and use those for mean rising limb (but use single peak events)
  #identify difference between peaks
  peak_dif_vec = rep(NA, length(peak_vec))
  peak_dif_vec[1] = 0
  for(i in c(2:length(peak_vec))){
    peak_dif_vec[i] = peak_vec[i]-peak_vec[i-1]
  }
  temp = data.frame(peak_vec, peak_dif_vec, Qobsdf[peak_vec,])
  #order dataframe by Q magnitude
  temp = temp[order(temp$Qobs, decreasing = T),]
  #select 5 highest events that have a peak rising time over longer than the 75% percentile of peak rising time (min 5 days)
  temp5 = temp[temp$peak_dif_vec>max(quantile(peak_dif_vec)[4], 5),][1:5,]
  #find valley points
  valley_vec = findValleys(Qobsdf[,2])-1 #is always off by one
  #calculate mean difference between 10 peaks and the valley right before
  mean_rising_time = mean(unlist(lapply(c(1:5), FUN = function(x){
    peakp = temp5$peak_vec[x]
    valleyp = tail(valley_vec[which(valley_vec<peakp)],1)
    peakp - valleyp
  })))
  #difference between peaks needs to be higher than mean rising time
  #highest peak is independent, two closest peaks need to be checked
  
  #do it iterativeley? select all POT above treshold, order by magnitude (highest first), then take that one and remove all peaks that are within the mean rising time
  #1. remove all peaks less than threshold
  temp1 = temp[temp$Qobs>POT_thresh,]
  #2. remove all peaks within peak rising time of a bigger peak
  #while loop (while tests if i is still bigger than the dimensions of the data frame)
  POTdf = temp1
  i = 1
  while(i<=dim(POTdf)[1]){
    #create time vector to catch all dates in +- rising time window of biggest peak
    date_vec = seq(POTdf[i,3]-(3*mean_rising_time), POTdf[i,3]+(3*mean_rising_time), by ="days")
    #which events are too close to the max peak (which is always the first)
    all_ind = as.character(POTdf[,3]) %in% as.character(date_vec) #change to character, as otherwise does not work
    #which one is the maximum peak
    keep_ind = which.max(as.character(POTdf[,3]) %in% as.character(date_vec)) #change to character, as otherwise does not work
    all_ind[keep_ind] = FALSE #reverse the maximum to keep this event
    #remove all other dependent events
    POTdf = POTdf[!all_ind,]
    i = i+1
  }
  POTdf = POTdf[,3:4]
  #events are ordered by magnitude, to get x events per year, calculate number of years and cut df accordingly 
  year_num = length(unique(years(Qobsdf[,1])))-1 #CAMELS starts with hydrological year, so remove 1
  POTdfy = POTdf[1:min((year_num*eventsperyear), length(POTdf[,1])),] #if there are not enough events, limit is set to maximum number of events
  return(POTdfy)
}



#normalise numeric variables
norm_func = function(x){
  (x-min(x, na.rm = T))/(max(x, na.rm = T)-min(x, na.rm = T))
}

#Coefficient of determination (R-squared)
r2func = function(obs, sim){cor(obs, sim)^2}

#change output from aggregate to a data frame
aggregatetodf = function(agg_out){cbind(agg_out[-ncol(agg_out)], agg_out[[ncol(agg_out)]])}


# Read in CAMELS catchment attributes -------------------------------------
camels_name = read.csv(file = paste0(CAMELS_path, "camels_attributes_v2.0/camels_name.txt"), sep = ";", stringsAsFactors = F, encoding = 'latin1')
camels_hydro = read.csv(file = paste0(CAMELS_path, "camels_attributes_v2.0/camels_hydro.txt"), sep = ";", stringsAsFactors = F, encoding = 'latin1')
camels_topo = read.csv(file =  paste0(CAMELS_path, "camels_attributes_v2.0/camels_topo.txt"), sep = ";", stringsAsFactors = F, encoding = 'latin1')
camels_clim = read.csv(file =  paste0(CAMELS_path, "camels_attributes_v2.0/camels_clim.txt"), sep = ";", stringsAsFactors = F, encoding = 'latin1')
camels_soil = read.csv(file =  paste0(CAMELS_path, "camels_attributes_v2.0/camels_soil.txt"), sep = ";", stringsAsFactors = F, encoding = 'latin1')
camels_vege = read.csv(file =  paste0(CAMELS_path, "camels_attributes_v2.0/camels_vege.txt"), sep = ";", stringsAsFactors = F, encoding = 'latin1')
camels_geol = read.csv(file =  paste0(CAMELS_path, "camels_attributes_v2.0/camels_geol.txt"), sep = ";", stringsAsFactors = F, encoding = 'latin1')
load(file = paste0(path_additional, "elongratio_CAMELS.Rdata"))
#from code file CAMELS_catchment_characteristics.R

camels_gaugeid = camels_name[,1]

#all continous variables
attr_cont = data.frame(camels_topo[,4:6], elongratio_CAMELS, camels_clim[,-c(1,9, 12)], camels_soil[,-1], camels_vege[,c(2:6)], camels_geol[,c(6:8)]) #removed log area log_area_topo
#normalise variables
attr_norm = data.frame(apply(attr_cont, 2, norm_func))


# Climate classification ---------------------------------------------------

clim_index = rep(NA, nrow(camels_clim))
clim_index[which(camels_clim$aridity<1)] = "Wet"
clim_index[which(camels_clim$aridity>=1)] = "Dry"
clim_index[which(camels_clim$frac_snow>=0.2)] = "Snow"

clim_indexdf = data.frame(camels_topo, clim_index, camels_clim)
clim_label_vec = c("Wet", "Dry", "Snow")


# Prepare flood classification ---------------------------------------------------

#Load available water storage data calculate from gNATSGO
AWC_CAMELS <- read.csv(paste0(path_additional, "aws0200NATSGO.csv"))
AWC_CAMELS = AWC_CAMELS[,c("HRU_ID", "MEAN")]
colnames(AWC_CAMELS) = c("gauge_id", "MEAN_AWS")
AWC_CAMELS[,"MEAN_AWS"] = AWC_CAMELS[, "MEAN_AWS"]*10 #unit is in cm, multiply by 10 for mm

#Model ouput files have climate and observed runoff all in  mm/d
model_output_qpath = list.dirs(paste0(CAMELS_path, "basin_timeseries_v1p2_modelOutput_daymet/model_output_daymet/model_output/flow_timeseries/daymet"))
qfiles =unlist(lapply(model_output_qpath, list.files, full.name = T, pattern = "05_model_output"))

# Flood classification ----------------------------------------------------
#Determine flood generating processes for each CAMELS catchment
event_mech_lists_POT3py = lapply(camels_gaugeid, function(catid){
  #read in discharge data
  qobstemp = read_qobs(catid, qfiles)
  #calculate Peaks-over-threshold
  POT_thresh = mean(qobstemp[,2])
  POT_temp = POT_func(qobstemp, POT_thresh, eventsperyear = 3)
  #read in climate data
  climate_in = read_climate(catid, qfiles)
  
  #calculate soil moisture and snowmelt
  #Tcrit changed to zero, because Addor, 2017 uses 0 degree to fraction snow calculation
  class_input_df_CAMELS = soil_snow_func(AWC_in = AWC_CAMELS[grep(catid, AWC_CAMELS[,1]),2], ET_mat = climate_in[,'ET'], rain_mat=climate_in[,'RAIM'], Temp_mat = climate_in[,'TAIR'], fdd_in = 2, Tcrit_in = 1, qobstemp[,1])
  #prepare flood peak data
  flood_df = list(list(catid, data.frame(AMAX_date = POT_temp[,1], AMAX = POT_temp[,2])))
  
  #create decision data frame
  event_df_list_CAMELS = event_classification_df(catid, day_thresh = 7, quant_thresh_rain = 0.9, quant_thresh_snow = 0.9, df_list = class_input_df_CAMELS, full_cat_list = catid, flood_df_in = flood_df, date_vec = qobstemp[,1])
  #classify process
  event_mech_list_out =event_classification_quantile(event_df_list_CAMELS, parts_rainsnow = 1/3, parts_fracextreme = 2/3, sat_thresh = 0.9)
})

save(event_mech_lists_POT3py, file = paste0(save_path, "event_mech_lists_POT3py.Rdata"))

event_mech_lists_POT1py = lapply(camels_gaugeid, function(catid){
  #read in discharge data
  qobstemp = read_qobs(catid, qfiles)
  #calculate Peaks-over-threshold
  POT_thresh = mean(qobstemp[,2])
  POT_temp = POT_func(qobstemp, POT_thresh, eventsperyear = 1)
  #read in climate data
  climate_in = read_climate(catid, qfiles)
  
  #calculate soil moisture and snowmelt
  #Tcrit changed to zero, because Addor, 2017 uses 0 degree to fraction snow calculation
  class_input_df_CAMELS = soil_snow_func(AWC_in = AWC_CAMELS[grep(catid, AWC_CAMELS[,1]),2], ET_mat = climate_in[,'ET'], rain_mat=climate_in[,'RAIM'], Temp_mat = climate_in[,'TAIR'], fdd_in = 2, Tcrit_in = 1, qobstemp[,1])
  #prepar flood peak data
  flood_df = list(list(catid, data.frame(AMAX_date = POT_temp[,1], AMAX = POT_temp[,2])))
  #create decision data frame
  event_df_list_CAMELS = event_classification_df(catid, day_thresh = 7, quant_thresh_rain = 0.9, quant_thresh_snow = 0.9, df_list = class_input_df_CAMELS, full_cat_list = catid, flood_df_in = flood_df, date_vec = qobstemp[,1])
  #classify process
  event_mech_list_out =event_classification_quantile(event_df_list_CAMELS, parts_rainsnow = 1/3, parts_fracextreme = 2/3, sat_thresh = 0.9)
})

save(event_mech_lists_POT1py, file = paste0(save_path, "event_mech_lists_POT1py.Rdata"))



# Summarise flood classification outcome ----------------------------------

#for three flood events per year
load(file = paste0(save_path, "event_mech_lists_POT3py.Rdata"))

level_vec = c("missingData", "other", "soilsat",  "rainfall", "longrainfall", "snowmelt", "rainandsnow")
level_vec_class= c("soilsat",  "rainfall", "longrainfall", "snowmelt", "rainandsnow")
count_mech_list = lapply(event_mech_lists_POT3py, function(event_mech){
  #Calculate percentage contribution, excluding missing and other events
  temp = event_mech[event_mech[,'mech_out'] %in% level_vec_class,]
  vec = factor(as.factor(temp[, 'mech_out']), levels = level_vec_class)
  table(vec)/length(vec)
})
count_mechdf = data.frame(do.call(rbind, count_mech_list))
#replace all NAs (due to missing data or class other) with zero to indicate zero process contribution
count_mechdf[is.na(count_mechdf)] = 0


#for one flood events per year
load(file = paste0(save_path, "event_mech_lists_POT1py.Rdata"))

count_mech_list = lapply(event_mech_lists_POT1py, function(event_mech){
  #Calculate percentage contribution, excluding missing and other events
  temp = event_mech[event_mech[,'mech_out'] %in% level_vec_class,]
  vec = factor(as.factor(temp[, 'mech_out']), levels = level_vec_class)
  table(vec)/length(vec)
})
count_mechdf1py = data.frame(do.call(rbind, count_mech_list))
count_mechdf1py[is.na(count_mechdf1py)] = 0


# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Probability distributions -----------------------------------------------
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

load(file = paste0(save_path, "event_mech_lists_POT3py.Rdata"))
load(file = paste0(save_path, "event_mech_lists_POT1py.Rdata"))

#Create dataframe where every event has all catchment attributes associated with it (POT1 and POT3)
event_based_magn = lapply(list(event_mech_lists_POT1py, event_mech_lists_POT3py), function(event_df){
  event_based_df_list = lapply(c(1:length(event_df)), function(cat){
    tempdf = event_df[[cat]]
    tempdf = tempdf[,c("mech_out", "sat_start", "Pmax", "P7", "S7", "flood_magn")]
    tempchar = attr_norm[cat,]
	#repeat attribute the number of events for that catchment
    tempchardf = tempchar[rep(seq_len(nrow(tempchar)),nrow(tempdf)),]
    tempclim = rep(clim_index[cat], nrow(tempdf))
    outdf = data.frame(tempdf, tempchardf, clim_class = tempclim, stringsAsFactors = T)
    return(outdf)
  })
  event_based_df = do.call(rbind, event_based_df_list)
  event_based_df = event_based_df[-which(event_based_df$mech_out == "noclass"),]
})

#count table
count_event_table = lapply(event_based_magn, function(df){
  aggregate(df$mech_out, by = list(df$clim_class), function(mech){
    vec = factor(as.factor(mech), levels = level_vec[-c(1,2)])
    table(vec)
  })
})
count_event_table
colSums(count_event_table[[2]][,-1])

count_event_df3 = aggregatetodf(count_event_table[[2]])
count_event_df1 = aggregatetodf(count_event_table[[1]])


#Calculate difference in distributions
magn_diff_list = lapply(event_based_magn, function(df){
  clim_plotdistlist_vsimple = lapply(colnames(attr_cont), function(tempvar){
    plotdist_list_clim = lapply(clim_label_vec, function(clim){
      clim_tempdf = df[df$clim_class == clim,]
      
	  #Process specific distribution
      plotdist_list = lapply(level_vec_class, function(process){
        clim_tempdf[which(clim_tempdf[,'mech_out']==process),tempvar]
      })
      #Process not specific distribution
	  plotdist_list_All = clim_tempdf[,tempvar]
	  #Calculate empirical distribution functions
      diff_result = do.call(rbind, lapply(plotdist_list[1:5], function(dist){
        ecdf1 = ecdf(dist)
        ecdf2 = ecdf(attr_norm[clim_index==clim,tempvar])
        out_mean= mean(ecdf2(seq(0,1, 0.01))-ecdf1(seq(0,1, 0.01)))
        out_abs = mean(abs(ecdf2(seq(0,1, 0.01))-ecdf1(seq(0,1, 0.01))))
        return(c(out_mean, out_abs))
      }))
      return(diff_result)
    })
    df_mean = do.call(rbind, lapply(plotdist_list_clim, function(x){x[,1]}))
    df_abs = do.call(rbind, lapply(plotdist_list_clim, function(x){x[,2]}))
    colnames(df_mean) = level_vec_class
    rownames(df_mean) = clim_label_vec
    colnames(df_abs) = level_vec_class
    rownames(df_abs) = clim_label_vec
    meltdf_mean = melt(df_mean)
    meltdf_abs = melt(df_abs)
    meltdf = data.frame(meltdf_mean, abs = meltdf_abs$value, attr = rep(tempvar, nrow(meltdf_mean)))
    return(meltdf)
  })
  return(clim_plotdistlist_vsimple)
})

char_group = c(rep("Topo", 4), rep("Clim", 9), rep("Soils", 11), rep("Vege", 5), rep("Geol", 3))
diff_df_3py = do.call(rbind, magn_diff_list[[2]])
colnames(diff_df_3py) = c("clim", "process", "mean_val", "abs_val", "attr")
#for plotting purposes define factor levels
diff_df_3py[,"attr"] = factor(diff_df_3py$attr, levels = rev(unique(diff_df_3py$attr)))
diff_df_3py[,"clim"] = factor(diff_df_3py$clim, levels = clim_label_vec)
diff_df_3py[,"process"] = factor(diff_df_3py$process, levels = level_vec_class)
diff_df_3py[,"char_group"] = factor(rep(char_group, each = 15), levels = unique(char_group))

count_event_df3_melt = melt(count_event_df3[c(2,3,1),], id = "Group.1") #reorder columns into wet, dry, snow
event_count = count_event_df3_melt[rep(1:nrow(count_event_df3_melt), length(unique(diff_df_3py$attr))),3]
diff_df_3py[,"event_count"] = event_count


#calculate for epy1
diff_df_1py = do.call(rbind, magn_diff_list[[1]])
colnames(diff_df_1py) = c("clim", "process", "mean_val", "abs_val", "attr")

diff_df_1py[,"attr"] = factor(diff_df_1py$attr, levels = rev(unique(diff_df_1py$attr)))
diff_df_1py[,"clim"] = factor(diff_df_1py$clim, levels = clim_label_vec)
diff_df_1py[,"process"] = factor(diff_df_1py$process, levels = level_vec_class)
diff_df_1py[,"char_group"] = factor(rep(char_group, each = 15), levels = unique(char_group))


count_event_df1_melt = melt(count_event_df1[c(2,3,1),], id = "Group.1") #reorder columns into wet, dry, snow
event_count = count_event_df1_melt[rep(1:nrow(count_event_df1_melt), length(unique(diff_df_1py$attr))),3]
diff_df_1py[,"event_count"] = event_count

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Random forest and interpretable machine learning ------------------------
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# RF cross-validation prediction all variables ----------------------------

set.seed(1) #for cross validation

clim_pred = lapply(clim_label_vec, function(clim){
  process_pred = lapply(level_vec_class, function(process_ind){
    tempdf = data.frame(process = count_mechdf[,process_ind], attr_cont[,-31]) #remove geol_porosity (includes NAs)
    tempdf = tempdf[clim_index==clim,]
    #draw cross validation samples
	fold_ind_clim = sample(rep(1:10, length.out = nrow(tempdf)), size = nrow(tempdf), replace = F)
    gaugeid_clim = camels_gaugeid[clim_index==clim]
	#iterate through cross-validation indices
    CV_rf <- lapply(1:10, function(fold){ 
      #fit model for 9/10 of data
	  fit <- randomForest(process~., data=tempdf[fold_ind_clim!=fold,])
	  #predict for remaining 1/10 of data
      predictions <- predict(fit, tempdf[fold_ind_clim==fold,], type = 'response')
      return(data.frame(Obs = tempdf$process[fold_ind_clim == fold], Pred = predictions, gauge_id = gaugeid_clim[fold_ind_clim==fold]))
    })
    CV_rf  = do.call(rbind, CV_rf)
    CV_rf = data.frame(CV_rf, process = rep(process_ind, nrow(CV_rf)))
    return(CV_rf)
  })
  process_pred_df = do.call(rbind, process_pred)
  process_pred_df = data.frame(process_pred_df, clim = rep(clim, nrow(process_pred_df)))
  return(process_pred_df)
})
clim_pred_df = do.call(rbind, clim_pred)


# Map differences in prediction -------------------------------------------

clim_pred_df_merge = merge(clim_pred_df, camels_topo[, c("gauge_id", "gauge_lat", "gauge_lon")])
clim_pred_df_merge[,"pred_diff"] = (clim_pred_df_merge$Pred-clim_pred_df_merge$Obs)^2
clim_pred_df_merge[,"pred_diff_norm"] = (clim_pred_df_merge$pred_diff/clim_pred_df_merge$Obs)
clim_pred_df_merge$pred_diff_norm[is.infinite(clim_pred_df_merge$pred_diff_norm)] = NA

# Accumulated local effects -------------------------------------------

char_group = c(rep("Topo", 4), rep("Clim", 9), rep("Soils", 11), rep("Vege", 5), rep("Geol", 2))#adjust length for removal of log area and geol_porostiy

#aggregate over all climates/processes
ALEdfagg = do.call(rbind, lapply(clim_label_vec, function(clim_in){
  out = do.call(rbind, lapply(level_vec_class, function(process_in){
    tempdf = data.frame(process = count_mechdf[,process_in], attr_cont[,-31])
    tempdf = tempdf[clim_index==clim_in,]
    #create random forest model
    rf = randomForest(process ~ ., data = tempdf)
    X =  tempdf[which(names(tempdf) != "process")]
    #Predict values
    model = Predictor$new(rf, data = X, y = tempdf$process)
    
    ALEdf = do.call(rbind, lapply(c(1:ncol(X)), function(i){
      #Calculate accumulated local effects
      ALEtemp = FeatureEffect$new(model, method = 'ale', feature = colnames(X)[i])
      ALEtempdf = ALEtemp$results
      data.frame(fval = ALEtempdf[,1], attr = rep(colnames(X)[i], nrow(ALEtempdf)))
    }))
    #take mean absolute value of accumulated effects for each attribute
    tempagg = data.frame(aggregatetodf(aggregate(ALEdf$fval, by = list(ALEdf$attr), function(x){mean(abs(as.numeric(x)))})), char_group, rep(process_in, length(char_group)), rep(clim_in, length(char_group)))
    colnames(tempagg) = c('attr','val', 'group', 'process', 'clim')
    
    #skewed distribution of input calculates less local effects, less reliable 
    temp_grid_no = do.call(rbind, lapply(unique(ALEdf$attr), function(attr){
      data.frame(attr = attr, ALE_nrow = nrow(ALEdf[ALEdf$attr== attr,]))
    }))
    tempagg = merge(tempagg, temp_grid_no, sort = F)
    tempagg[,'val'] = norm_func(tempagg[,'val'])
    return(tempagg)
  }))
}))

ALEdfagg[,"group"] = factor(ALEdfagg$group, levels = unique(char_group))
head(ALEdfagg)


# ALE for POT1
#aggregate over all climates/processes
ALEdfagg_POT1 = do.call(rbind, lapply(clim_label_vec, function(clim_in){
  out = do.call(rbind, lapply(level_vec_class, function(process_in){
    tempdf = data.frame(process = count_mechdf1py[,process_in], attr_cont[,-31])
    tempdf = tempdf[clim_index==clim_in,]
    #create random forest model
    rf = randomForest(process ~ ., data = tempdf)
    X =  tempdf[which(names(tempdf) != "process")]
    #Predict values
    model = Predictor$new(rf, data = X, y = tempdf$process)
    
    ALEdf = do.call(rbind, lapply(c(1:ncol(X)), function(i){
      #Calculate accumulated local effects
      ALEtemp = FeatureEffect$new(model, method = 'ale', feature = colnames(X)[i])
      ALEtempdf = ALEtemp$results
      data.frame(fval = ALEtempdf[,1], attr = rep(colnames(X)[i], nrow(ALEtempdf)))
    }))
    #take mean absolute value of accumulated effects for each attribute
    tempagg = data.frame(aggregatetodf(aggregate(ALEdf$fval, by = list(ALEdf$attr), function(x){mean(abs(as.numeric(x)))})), char_group, rep(process_in, length(char_group)), rep(clim_in, length(char_group)))
    colnames(tempagg) = c('attr','val', 'group', 'process', 'clim')
    #skewed distribution of input calculates less local effects, less reliable 
    temp_grid_no = do.call(rbind, lapply(unique(ALEdf$attr), function(attr){
      data.frame(attr = attr, ALE_nrow = nrow(ALEdf[ALEdf$attr== attr,]))
    }))
    tempagg = merge(tempagg, temp_grid_no, sort = F)
    tempagg[,'val'] = norm_func(tempagg[,'val'])
    return(tempagg)
  }))
}))

ALEdfagg_POT1[,"group"] = factor(ALEdfagg_POT1$group, levels = unique(char_group))



# Ablation analysis  -------------------------------------------



attr_RF =attr_cont[, -31] #remove log area and geol_porosity (includes NAs)
char_group = c(rep("Topo", 4), rep("Clim", 9), rep("Soils", 11), rep("Vege", 5), rep("Geol", 2))

comb_list = lapply(c(1:5), function(i){
  combinations(n = 5, r = i, v = unique(char_group), repeats.allowed = F)
})
#make string vector for all combinations
comb_vec = unlist(lapply(comb_list, function(x){apply(x,1,paste,collapse=" ")}))

#make list of attribute data frames according to combinations of attributes
attr_RF_list = do.call(c,lapply(comb_list, function(combmat){
  apply(combmat, 1, function(x){attr_RF[,(char_group %in% x)]})
}))

#
#Iterate trough all possible combinations of attribute groups (with removal)
s_time = Sys.time()
# Calculate the number of cores
no_cores <- detectCores() - 3
# Initiate cluster
cl <- makeCluster(no_cores)
clusterExport(cl, varlist = c("attr_RF_list", "count_mechdf", "level_vec_class", "clim_index",  "randomForest", "clim_label_vec"))
RF_diff_attr = parLapply(cl, attr_RF_list, function(attr_df){
  set.seed(1)
  clim_pred = lapply(clim_label_vec, function(clim){
    process_pred = lapply(level_vec_class, function(process_ind){
      tempdf = data.frame(process = count_mechdf[,process_ind], attr_df)
      tempdf = tempdf[clim_index==clim,]
      fold_ind_clim = sample(rep(1:10, length.out = nrow(tempdf)), size = nrow(tempdf), replace = F)
      CV_rf <- lapply(1:10, function(fold){ 
        fit <- randomForest::randomForest(process~., data=tempdf[fold_ind_clim!=fold,])
        predictions <- predict(fit, tempdf[fold_ind_clim==fold,], type = 'response')
        return(data.frame(Obs = tempdf$process[fold_ind_clim == fold], Pred = predictions))
      })
      CV_rf  = do.call(rbind, CV_rf)
      CV_rf = data.frame(CV_rf, process = rep(process_ind, nrow(CV_rf)))
      return(CV_rf)
    })
    process_pred_df = do.call(rbind, process_pred)
    process_pred_df = data.frame(process_pred_df, clim = rep(clim, nrow(process_pred_df)))
    return(process_pred_df)
  })
  clim_pred_df = do.call(rbind, clim_pred)
  return(clim_pred_df)
})
stopCluster(cl)
e_time = Sys.time()
e_time-s_time
save(RF_diff_attr, file = paste0(save_path, "RF_diff_attr.Rdata"))


load(file = paste0(save_path, "RF_diff_attr.Rdata"))
#Evaluate random forests through KGE and r2
eval_diff_attr = lapply(RF_diff_attr, function(clim_pred_df){
  lapply(level_vec_class, function(process){
    tempdf = clim_pred_df[clim_pred_df$process == process,]
    out = unlist(lapply(clim_label_vec, function(clim){
      tempdf = tempdf[tempdf$clim == clim,]
      KGE_out = unlist(KGE(tempdf[,2], tempdf[,1], out.type = "full"))
      round(c(r2_out = r2func(tempdf[,2], tempdf[,1]), KGE_out),4)
    }))
    return(out)
  })
})

r2_eval = lapply(eval_diff_attr, function(eval_list){
  do.call(rbind, lapply(eval_list, function(x){x[grepl("^r2_out", names(x))]})) 
})

#function to label which group got removed
remove_group_func = function(bestdf, worstdf){
  a = strsplit(as.character(bestdf[1,2]), " ")[[1]]
  b = strsplit(as.character(worstdf[1,2]), " ")[[1]]
  remove_group = a[!(a %in% b)]
  return(remove_group)
}
#Remove group contributing most (='best')
comb_length = unlist(lapply(strsplit(comb_vec, " "), length))
ablation_best_full = do.call(rbind, lapply(c(1:5), function(process_in){
  process = level_vec_class[process_in]
  do.call(rbind, lapply(c(1:3), function(clim_in){
    clim = clim_label_vec[clim_in]
    #select R2 for process, clim combination
    r2_eval_temp = unlist(lapply(r2_eval, function(x)x[process_in,clim_in]))
    #how many groups are removed
    r2comb = data.frame(r2_eval_temp, comb_vec, leng = 6-comb_length)
    out = list()
    for(i in c(1:5)){
      #R2 for combination with most groups
      allr2 = r2comb[nrow(r2comb),]
      #which one is the worst combination (= max difference) with one less group as part of the combination
      worst_next = r2comb[r2comb$leng == i+1,][which.max(allr2[1,1]-r2comb[r2comb$leng == i+1,1]),]
      #which group got removed
      remove_group = remove_group_func(allr2, worst_next)
      #new starting df
      r2comb = r2comb[!grepl(remove_group, r2comb[,2]),]
      out[[i]] = data.frame(allr2, remove_group)
    }
    ablation_best = do.call(rbind, out)
    #instead of remove group, list removed group
    ablation_best$remove_group = c('All', as.character(ablation_best$remove_group[-5]))
    ablation_best[,"comb_vec"] = factor(ablation_best$comb_vec, levels = ablation_best$comb_vec)
    ablation_best = data.frame(ablation_best, process = rep(process, nrow(ablation_best)), clim = rep(clim, nrow(ablation_best)))
    return(ablation_best)
  }))
}))


#Remove group contributing least (='worst')
ablation_rev_full = do.call(rbind, lapply(c(1:5), function(process_in){
  process = level_vec_class[process_in]
  do.call(rbind, lapply(c(1:3), function(clim_in){
    clim = clim_label_vec[clim_in]
    r2_eval_temp = unlist(lapply(r2_eval, function(x)x[process_in,clim_in]))
    r2comb = data.frame(r2_eval_temp, comb_vec, leng = 6-comb_length)
    out = list()
    for(i in c(1:5)){
      allr2 = r2comb[nrow(r2comb),]
      worst_next = r2comb[r2comb$leng == i+1,][which.min(allr2[1,1]-r2comb[r2comb$leng == i+1,1]),]
      remove_group = remove_group_func(allr2, worst_next)
      r2comb = r2comb[!grepl(remove_group, r2comb[,2]),]
      out[[i]] = data.frame(allr2, remove_group)
    }
    ablation_rev = do.call(rbind, out)
    #instead of remove group, list removed group
    ablation_rev$remove_group = c('All', as.character(ablation_rev$remove_group[-5]))
    ablation_rev[,"comb_vec"] = factor(ablation_rev$comb_vec, levels = ablation_rev$comb_vec)
    ablation_rev = data.frame(ablation_rev, process = rep(process, nrow(ablation_rev)), clim = rep(clim, nrow(ablation_rev)))
    return(ablation_rev)
  }))
}))
ablation_best_full[,"remove_group"] = factor(ablation_best_full$remove_group, levels = c('All', unique(char_group)))
ablation_rev_full[,"remove_group"] = factor(ablation_rev_full$remove_group, levels = c('All', unique(char_group)))


# repeat for POT1 -------------------
#Iterate trough all possible combinations of attribute groups (with removal)
s_time = Sys.time()
# Calculate the number of cores
no_cores <- detectCores() - 3
# Initiate cluster
cl <- makeCluster(no_cores)
clusterExport(cl, varlist = c("attr_RF_list", "count_mechdf1py", "level_vec_class", "clim_index",  "randomForest", "clim_label_vec"))
RF_diff_attr_POT1 = parLapply(cl, attr_RF_list, function(attr_df){
  set.seed(1)
  clim_pred = lapply(clim_label_vec, function(clim){
    process_pred = lapply(level_vec_class, function(process_ind){
      tempdf = data.frame(process = count_mechdf1py[,process_ind], attr_df)
      tempdf = tempdf[clim_index==clim,]
      fold_ind_clim = sample(rep(1:10, length.out = nrow(tempdf)), size = nrow(tempdf), replace = F)
      CV_rf <- lapply(1:10, function(fold){ 
        fit <- randomForest::randomForest(process~., data=tempdf[fold_ind_clim!=fold,])
        predictions <- predict(fit, tempdf[fold_ind_clim==fold,], type = 'response')
        return(data.frame(Obs = tempdf$process[fold_ind_clim == fold], Pred = predictions))
      })
      CV_rf  = do.call(rbind, CV_rf)
      CV_rf = data.frame(CV_rf, process = rep(process_ind, nrow(CV_rf)))
      return(CV_rf)
    })
    process_pred_df = do.call(rbind, process_pred)
    process_pred_df = data.frame(process_pred_df, clim = rep(clim, nrow(process_pred_df)))
    return(process_pred_df)
  })
  clim_pred_df = do.call(rbind, clim_pred)
  return(clim_pred_df)
})
stopCluster(cl)
e_time = Sys.time()
e_time-s_time
save(RF_diff_attr_POT1, file = paste0(save_path, "RF_diff_attr_POT1.Rdata"))



load(file = paste0(save_path, "RF_diff_attr_POT1.Rdata"))
#Evaluate random forests through KGE and r2
eval_diff_attr_POT1 = lapply(RF_diff_attr_POT1, function(clim_pred_df){
  lapply(level_vec_class, function(process){
    tempdf = clim_pred_df[clim_pred_df$process == process,]
    out = unlist(lapply(clim_label_vec, function(clim){
      tempdf = tempdf[tempdf$clim == clim,]
      KGE_out = unlist(KGE(tempdf[,2], tempdf[,1], out.type = "full"))
      round(c(r2_out = r2func(tempdf[,2], tempdf[,1]), KGE_out),4)
    }))
    return(out)
  })
})

r2_eval_POT1 = lapply(eval_diff_attr_POT1, function(eval_list){
  do.call(rbind, lapply(eval_list, function(x){x[grepl("^r2_out", names(x))]})) 
})


#Remove group contributing most (='best')
comb_length = unlist(lapply(strsplit(comb_vec, " "), length))
ablation_best_full_POT1 = do.call(rbind, lapply(c(1:5), function(process_in){
  process = level_vec_class[process_in]
  do.call(rbind, lapply(c(1:3), function(clim_in){
    clim = clim_label_vec[clim_in]
    #select R2 for process, clim combination
    r2_eval_temp = unlist(lapply(r2_eval_POT1, function(x)x[process_in,clim_in]))
    #how many groups are removed
    r2comb = data.frame(r2_eval_temp, comb_vec, leng = 6-comb_length)
    out = list()
    for(i in c(1:5)){
      #R2 for combination with most groups
      allr2 = r2comb[nrow(r2comb),]
      #which one is the worst combination (= max difference) with one less group as part of the combination
      worst_next = r2comb[r2comb$leng == i+1,][which.max(allr2[1,1]-r2comb[r2comb$leng == i+1,1]),]
      #which group got removed
      remove_group = remove_group_func(allr2, worst_next)
      #new starting df
      r2comb = r2comb[!grepl(remove_group, r2comb[,2]),]
      out[[i]] = data.frame(allr2, remove_group)
    }
    ablation_best = do.call(rbind, out)
    #instead of remove group, list removed group
    ablation_best$remove_group = c('All', as.character(ablation_best$remove_group[-5]))
    ablation_best[,"comb_vec"] = factor(ablation_best$comb_vec, levels = ablation_best$comb_vec)
    ablation_best = data.frame(ablation_best, process = rep(process, nrow(ablation_best)), clim = rep(clim, nrow(ablation_best)))
    return(ablation_best)
  }))
}))


#Remove group contributing least (='worst')
ablation_rev_full_POT1 = do.call(rbind, lapply(c(1:5), function(process_in){
  process = level_vec_class[process_in]
  do.call(rbind, lapply(c(1:3), function(clim_in){
    clim = clim_label_vec[clim_in]
    r2_eval_temp = unlist(lapply(r2_eval_POT1, function(x)x[process_in,clim_in]))
    r2comb = data.frame(r2_eval_temp, comb_vec, leng = 6-comb_length)
    out = list()
    for(i in c(1:5)){
      allr2 = r2comb[nrow(r2comb),]
      worst_next = r2comb[r2comb$leng == i+1,][which.min(allr2[1,1]-r2comb[r2comb$leng == i+1,1]),]
      remove_group = remove_group_func(allr2, worst_next)
      r2comb = r2comb[!grepl(remove_group, r2comb[,2]),]
      out[[i]] = data.frame(allr2, remove_group)
    }
    ablation_rev = do.call(rbind, out)
    #instead of remove group, list removed group
    ablation_rev$remove_group = c('All', as.character(ablation_rev$remove_group[-5]))
    ablation_rev[,"comb_vec"] = factor(ablation_rev$comb_vec, levels = ablation_rev$comb_vec)
    ablation_rev = data.frame(ablation_rev, process = rep(process, nrow(ablation_rev)), clim = rep(clim, nrow(ablation_rev)))
    return(ablation_rev)
  }))
}))
ablation_best_full_POT1[,"remove_group"] = factor(ablation_best_full_POT1$remove_group, levels = c('All', unique(char_group)))
ablation_rev_full_POT1[,"remove_group"] = factor(ablation_rev_full_POT1$remove_group, levels = c('All', unique(char_group)))




#-------------------





# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Visualisation setup -----------------------------------------------------
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#Define map outline and design
USmapdata = map_data("state")
USmap = ggplot() + geom_polygon(data = USmapdata, aes(x = long, y = lat, group = group), colour = "grey70", fill = "white", size = 0.1)
background = theme(panel.grid.major = element_line(colour = "grey90"), panel.background = element_rect(fill = "white", colour = 'transparent'),legend.text = element_text(size = 11), legend.title = element_text(size = 12), legend.key=element_blank(), legend.position = "bottom", legend.background = element_rect(fill="white"), axis.line=element_blank(), axis.title=element_blank(), strip.background = element_rect(fill  = 'white'), strip.text = element_text(size = 11))#+ coord_fixed(ratio = 2) #+ 
#Color scale for flood generating processes
#Colorblind safe
# taken from: https://personal.sron.nl/~pault/#sec:qualitative
Excessrain_col = '#228833'
Rainfall_col = '#EE6677'
Longrainfall_col = '#CCBB44'
Rainsnow_col = '#66CCEE'
Snowmelt_col = '#AA3377'
Other_col = '#BBBBBB'
Missing_col = 'white'
col_info = scale_colour_manual(name = "Flood process", labels = c("Excess rainfall", "Short rainfall",  "Long rainfall", "Snowmelt","Rainfall/Snowmelt", "All" ), values = c(Excessrain_col, Rainfall_col, Longrainfall_col, Snowmelt_col,Rainsnow_col, 'black'))

#facet strip design
strip_theme = theme(strip.background = element_blank(), strip.text.y = element_text(angle = 0), strip.text = element_text(face = 'bold'))


#Color scale for climate classes
clim_col_vec = c("#E7B800","#FC4E07", "#00AFBB")


#labels suitable for ggplot facets 
process_fulllabel = c('soilsat'= "Excess rainfall", 'rainfall' = "Short rainfall",  'longrainfall' = "Long rainfall", 'snowmelt' = "Snowmelt",'rainandsnow' = "Rainfall/Snowmelt")


#CAMELS attribute labels
CAMELS_labels = c(strsplit(readLines('C:/Users/ls16959/Data/streamflow_countries/US/Camels/CAMELS_parameter_labels.txt', warn = F), split = ', ')[[1]])

attr_fulllabel = CAMELS_labels
names(attr_fulllabel) = colnames(attr_cont)


#Shift legend to empty facet function
#Source: https://stackoverflow.com/questions/54438495/shift-legend-into-empty-facets-of-a-faceted-plot-in-ggplot2
shift_legend2 <- function(p) {
  # ...
  # to grob
  gp <- ggplotGrob(p)
  facet.panels <- grep("^panel", gp[["layout"]][["name"]])
  empty.facet.panels <- sapply(facet.panels, function(i) "zeroGrob" %in% class(gp[["grobs"]][[i]]))
  empty.facet.panels <- facet.panels[empty.facet.panels]
  
  # establish name of empty panels
  empty.facet.panels <- gp[["layout"]][empty.facet.panels, ]
  names <- empty.facet.panels$name
  # example of names:
  #[1] "panel-3-2" "panel-3-3"
  
  # now we just need a simple call to reposition the legend
  reposition_legend(p, 'center', panel=names)
}


# Plot Climate classification ---------------------------------------------------

p1 = ggplot(clim_indexdf, aes(aridity, frac_snow, col = clim_index))+geom_point(alpha = 0.6)+theme_minimal()+scale_colour_manual(limits = clim_label_vec,  values = clim_col_vec, name = 'Climate type')+xlab("Aridity")+ylab("Fraction of snow")+geom_hline(yintercept = 0.2, linetype = 2, col = "grey50")+geom_segment(aes(x = 1, y = 0, xend = 1, yend = 0.2), col = "grey50", linetype = 2)+theme(legend.position = c(0.8,0.8), panel.grid = element_blank())
p1
point_data = geom_point(data = clim_indexdf, aes(x=gauge_lon, y=gauge_lat, col = as.factor(clim_index)),na.rm = TRUE, size = 2, alpha = 0.6)
p2 = USmap+background+point_data+scale_colour_manual(limits = clim_label_vec,  values = clim_col_vec, name = '')+coord_map()+theme(legend.position = "none")
p2
ggsave(filename = paste0(plot_path, "CAMELS_simpleMAP.pdf"), dpi = 'retina', width = 10, height =7)

p1 + p2 +plot_layout(widths = c(1,2))+plot_annotation(tag_levels = 'A')

ggsave(filename = paste0(plot_path, "CAMELS_climatetype.pdf"), dpi = 'retina', width = 10, height =4)

# Plot gNATSGO available water storage ---------------------------------------------------
AWCplotdf = merge(AWC_CAMELS, camels_topo)
point_data = geom_point(data = AWCplotdf, aes(x=gauge_lon, y=gauge_lat, col = MEAN_AWS),na.rm = TRUE, size = 2)
USmap+background+point_data+coord_map()+scale_color_gradientn(colours=pals::brewer.blues(100), name = 'Availabel Water Storage')+guides(col = guide_colourbar(barwidth = 20))
ggsave(filename = paste0(plot_path, "CAMELS_AWS.pdf"), dpi = 'retina', width = 7, height =5)


# Camels attribute data overview plots -------------------------------------------------
#Data histograms
attr_melt = melt(attr_cont)
ggplot(attr_melt, aes(value))+geom_histogram()+xlab("Attribute value")+ylab("Count")+theme_bw()+theme(panel.grid = element_blank())+facet_wrap(~variable, scales = "free")
ggsave(filename = paste0(plot_path, "CAMELS_attr_hist.pdf"), dpi = 'retina', width = 10, height =7)

#Correlation of attributes with hydrological indices
attr_hyd = data.frame(attr_cont, camels_hydro[,-1])
corr_df = lapply(clim_label_vec, function(clim){
  tempdf = cor(attr_hyd[clim_index==clim,], use = "complete.obs", method = "spearman")
  tempdf = tempdf[1:32,33:45]
  tempmelt = melt(tempdf)
  out = data.frame(tempmelt, clim = rep(clim, nrow(tempmelt)))
  return(out)
})
corr_df = do.call(rbind, corr_df)
corr_df[,"attr1"] = factor(corr_df$X1, levels = rev(unique(corr_df$X1)))
corr_df[,"attr2"] = factor(corr_df$X2, levels = unique(corr_df$X2))

char_group = c(rep("Topo", 4), rep("Clim", 9), rep("Soils", 11), rep("Vege", 5), rep("Geol", 3))
corr_df[,"group"] = factor(rep(char_group, 3*length(unique(corr_df$attr2))), levels = unique(char_group))


ggplot(corr_df, aes(attr2, attr1, col = value, size = abs(value)))+geom_point()+scale_colour_gradient2(limits = c(-1, 1), name = "Cor", high = '#08306B')+theme_bw()+theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title = element_blank(), panel.border = element_rect(colour = 'grey80'))+strip_theme+scale_size(range = c(0.3, 4), name = "")+geom_vline(xintercept = c(6.5, 7.5))+facet_grid(group~clim, scales = "free", space = "free")

ggsave(filename = paste0(plot_path, "CAMELS_attr_hyd_corr_clim.pdf"), dpi = 'retina', width = 8, height = 6)


#Correlation of all attributes
tempdf = cor(attr_cont, use = "complete.obs", method = "spearman")
tempmelt = melt(tempdf)
tempmelt[,'X1_char_group'] = rep(char_group, times = 32)
tempmelt[,'X2_char_group'] = rep(char_group, each = 32)
tempmelt[,"attr1"] = factor(tempmelt$X1, levels = unique(tempmelt$X1))
tempmelt[,"attr2"] = factor(tempmelt$X2, levels = rev(unique(tempmelt$X2)))


ggplot(tempmelt, aes(attr1, attr2, col = value, size = abs(value)))+geom_point()+scale_colour_gradient2(limits = c(-1, 1), name = 'Cor')+theme_bw()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0), axis.title = element_blank())+scale_x_discrete(position = "top", labels = CAMELS_labels)+scale_y_discrete(labels = rev(CAMELS_labels))+scale_size(range = c(0.3, 4), name = '')+geom_hline(yintercept = c(3.5, 8.5, 19.5, 28.5))+geom_vline(xintercept = c(4.5, 13.5, 24.5, 29.5))
ggsave(filename = paste0(plot_path, "CAMELS_attr_corr_all_presentation.pdf"), dpi = 'retina', width = 7, height = 6) # width = 7, height = 6




#Correlation of all attributes split by climate
corr_df = lapply(clim_label_vec, function(clim){
  tempdf = cor(attr_cont[clim_index==clim,], use = "complete.obs", method = "spearman")
  tempmelt = melt(tempdf)
  out = data.frame(tempmelt, clim = rep(clim, nrow(tempmelt)))
  return(out)
})
corr_df = do.call(rbind, corr_df)
corr_df[,"attr1"] = factor(corr_df$X1, levels = unique(corr_df$X1))
corr_df[,"attr2"] = factor(corr_df$X2, levels = rev(unique(corr_df$X2)))


ggplot(corr_df, aes(attr1, attr2, col = value))+geom_point()+scale_colour_gradient2(limits = c(-1, 1), name = expression(rho))+theme_bw()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0), axis.title = element_blank())+scale_x_discrete(position = "top")+geom_hline(yintercept = c(3.5, 8.5, 19.5, 28.5))+geom_vline(xintercept = c(4.5, 13.5, 24.5, 29.5))+facet_wrap(~clim, nrow = 3, strip.position = "right")+strip_theme
ggsave(filename = paste0(plot_path, "CAMELS_attr_corr_clim.pdf"), dpi = 'retina', width = 6, height = 12)

#Alternative Plotting Option
attr_cont_temp = attr_cont[,c(ncol(attr_cont):1)]
dat <- data.frame(x = seq(attr_cont_temp), y = seq(attr_cont_temp))
#Wet
gcor = ggcorr(attr_cont_temp[which(clim_index == "Wet"),], low = muted("red"), mid = "white",high = muted("blue"), midpoint = 0, geom = "tile", hjust = 1, layout.exp = 13, color = NA, method = c("complete.obs", "spearman"),  name = expression(rho))
cor1 = gcor + geom_text(data=dat, aes(x, y, label=rev(CAMELS_labels)), nudge_x = 13, hjust=1, nudge_y = 0.3)+ ggplot2::labs(title = "Wet")
#dry
gcor = ggcorr(attr_cont_temp[which(clim_index == "Dry"),], low = muted("red"), mid = "white",high = muted("blue"), midpoint = 0, geom = "tile", hjust = 1, layout.exp = 13, color = NA, method = c("complete.obs", "spearman"),  name = expression(rho))
cor2 = gcor + geom_text(data=dat, aes(x, y, label=rev(CAMELS_labels)), nudge_x = 13, hjust=1, nudge_y = 0.3)+ ggplot2::labs(title = "Dry")
#Snow
gcor = ggcorr(attr_cont_temp[which(clim_index == "Snow"),], low = muted("red"), mid = "white",high = muted("blue"), midpoint = 0, geom = "tile", hjust = 1, layout.exp = 13, color = NA, method = c("complete.obs", "spearman"),  name = expression(rho))
cor3 = gcor + geom_text(data=dat, aes(x, y, label=rev(CAMELS_labels)), nudge_x = 13, hjust=1, nudge_y = 0.3)+ ggplot2::labs(title = "Snow")

cor1/cor2/cor3 +plot_layout(guides = "collect")

ggsave(filename = paste0(plot_path, "CAMELS_attr_corr_clim.pdf"), dpi = 'retina', width = 7, height = 14)



# Plot flood process distribution ---------------------------------------------------

process_mapdf = data.frame(camels_topo[,c("gauge_lat", "gauge_lon")], count_mechdf)
melt_process_mapdf = melt(process_mapdf, id = c("gauge_lat", "gauge_lon"))

#Process map
point_data = geom_point(data = melt_process_mapdf, aes(x=gauge_lon, y=gauge_lat, col = value*100),na.rm = TRUE, size = 2)
p =USmap+background+point_data+scale_color_gradientn(colours=pals::brewer.blues(100), name = 'Process contribution [%]', limits = c(0,100), breaks = c(0, 25, 50, 75, 100))+coord_map()+strip_theme+facet_wrap(~variable, labeller = labeller(variable = process_fulllabel), nrow = 2)+guides(col = guide_colourbar(barwidth = 12.5, title.position = 'top'))
#ggsave(filename = paste0(plot_path, "CAMELS_process_distribution.pdf"), dpi = 'retina', width = 10, height =5)
p1 = shift_legend2(p)
ggsave(filename = paste0(plot_path, "CAMELS_process_distribution.pdf"), plot = p1, dpi = 'retina', width = 10, height =4)

#Process barplot
tempmelt = do.call(rbind, lapply(clim_label_vec, function(clim){
  melt(data.frame(cat_no = c(1:nrow(count_mechdf[clim_index == clim,])), count_mechdf[clim_index == clim,], clim_index = clim_index[clim_index == clim]), id.vars = c('cat_no', 'clim_index'))
}))

p_contr = ggplot(tempmelt, aes(cat_no, 100*value, fill = variable))+geom_bar(stat = 'identity', width=1)+scale_fill_manual(name = "Flood process", labels = process_fulllabel, values = c(Excessrain_col, Rainfall_col, Longrainfall_col, Snowmelt_col,Rainsnow_col))+ylab('Process contribution [%]')+xlab('Catchments')+theme_minimal()+
  theme(panel.grid.minor  = element_blank(), panel.grid.major = element_blank())+
  strip_theme+facet_wrap(~clim_index, scales = 'free_x')
ggsave(filename = paste0(plot_path, "CAMELS_process_distribution_barplot.pdf"), dpi = 'retina', width = 10, height =4)


#Processs event numbers
count_event_table_POT3 = count_event_table[[2]][,-1]
rownames(count_event_table_POT3) = count_event_table[[2]][,1]
count_event_table_POT3_melt = melt(count_event_table_POT3)
count_event_table_POT3_melt[, "X1"] = factor(count_event_table_POT3_melt$X1, levels = clim_label_vec)
count_event_table_POT3_melt[, "X2"] = factor(count_event_table_POT3_melt$X2, levels = level_vec_class)



p_events = ggplot(count_event_table_POT3_melt, aes(X2, value, fill = X2))+geom_bar(stat = 'identity')+scale_fill_manual(name = "Flood process", labels = process_fulllabel, values = c(Excessrain_col, Rainfall_col, Longrainfall_col, Snowmelt_col,Rainsnow_col))+theme_bw()+theme(panel.border = element_blank(), panel.grid.minor  = element_blank(), panel.grid.major.x = element_blank(), axis.title.x = element_blank(), legend.position = 'none', axis.text.x = element_blank(), axis.ticks.x = element_blank())+strip_theme+facet_wrap(~X1,ncol = 3)+ylab('No. events')


p_contr/p_events + plot_layout(guides = "collect")+plot_annotation(tag_levels = 'A')

ggsave(filename = paste0(plot_path, "CAMELS_process_distribution_barplot.pdf"), dpi = 'retina', width = 10, height =6)




#for POT1
tempmelt = do.call(rbind, lapply(clim_label_vec, function(clim){
  melt(data.frame(cat_no = c(1:nrow(count_mechdf1py[clim_index == clim,])), count_mechdf1py[clim_index == clim,], clim_index = clim_index[clim_index == clim]), id.vars = c('cat_no', 'clim_index'))
}))

ggplot(tempmelt, aes(cat_no, 100*value, fill = variable))+geom_bar(stat = 'identity', width=1)+scale_fill_manual(name = "Flood process", labels = process_fulllabel, values = c(Excessrain_col, Rainfall_col, Longrainfall_col, Snowmelt_col,Rainsnow_col))+ylab('Process contribution [%]')+xlab('Catchments')+theme_minimal()+theme(legend.position = 'bottom')+strip_theme+facet_wrap(~clim_index, scales = 'free_x')
ggsave(filename = paste0(plot_path, "CAMELS_process_distribution_barplot_POT1.pdf"), dpi = 'retina', width = 10, height =4)


# Plot probability density distribution ---------------------------------------------------
strip_theme2 = theme(strip.background = element_blank(), strip.text = element_text(face = 'bold'), strip.text.y = element_text(angle = 0), panel.grid.minor = element_blank(), panel.grid.major = element_blank())

magn_plot_vsimple = lapply(event_based_magn, function(df){
  clim_plotdistlist_vsimple = lapply(colnames(attr_cont), function(tempvar){
    plotdist_list_clim = lapply(clim_label_vec, function(clim){
      clim_tempdf = df[df$clim_class == clim,]
      
      plotdist_list = lapply(level_vec_class, function(process){
        clim_tempdf[which(clim_tempdf[,'mech_out']==process),tempvar]
      })
      #use attribute distribution as 'All'
      plotdist_list[["All"]] = attr_norm[clim_index==clim,tempvar]
      level_vec_ext = c(level_vec_class, "All")
      length_vec = lapply(plotdist_list, length)
      df <- data.frame(x = unlist(plotdist_list), ggg=factor(rep(1:length(plotdist_list), length_vec)), clim_class = factor(rep(clim, length(unlist(plotdist_list)))), processname = factor(rep(level_vec_ext, length_vec)))
      return(df)
    })
    df = do.call(rbind,plotdist_list_clim)
    df = data.frame(df, attr = rep(tempvar, nrow(df)))
    return(df)
  })
})

xlab_info = xlab("Normalised attribute")
ylab_info = ylab("ECDF")
format_info = theme(legend.position = "bottom")

dist_list_1py = magn_plot_vsimple[[1]]
dist_list_3py = magn_plot_vsimple[[2]]

#POT1 and POT3 distribution side by side ----
pdf(paste0(plot_path, "CAMELS_dist_comparison_All.pdf"),width = 8, height =6)
for(i in c(1:length(dist_list_1py))){
  plotdf_temp = do.call(rbind, dist_list_1py[i])
  p1 = ggplot(plotdf_temp, aes(x, colour = ggg))+stat_ecdf(size = 0.8)+col_info+xlab_info+ylab_info+theme_bw()+scale_x_continuous(breaks = c(0.5, 1))+format_info+facet_wrap(~clim_class)+ggtitle(colnames(attr_cont)[i])
  
  plotdf_temp = do.call(rbind, dist_list_3py[i])
  p3 = ggplot(plotdf_temp, aes(x, colour = ggg))+stat_ecdf(size = 0.8)+col_info+xlab_info+ylab_info+theme_bw()+scale_x_continuous(breaks = c(0.5, 1))+format_info+facet_wrap(~clim_class)+ggtitle(colnames(attr_cont)[i])
  multiplot(p1, p3, cols = 1)
}
dev.off()


#write plots POT1py to files ----

#loco topo
plotdf_temp = do.call(rbind, dist_list_1py[c(1:4)])
ggplot(plotdf_temp, aes(x, colour = ggg))+stat_ecdf(size = 0.8)+col_info+xlab_info+ylab_info+theme_bw()+scale_x_continuous(breaks = c(0.5, 1))+format_info+strip_theme+facet_grid(attr~clim_class)
ggsave(filename = paste0(plot_path, "CAMELS_dist_comparison_locotopo1.pdf"), dpi = 'retina', width = 7, height =7)
#clim
plotdf_temp = do.call(rbind, dist_list_1py[c(5:13)])
ggplot(plotdf_temp, aes(x, colour = ggg))+stat_ecdf(size = 0.8)+col_info+xlab_info+ylab_info+theme_bw()+scale_x_continuous(breaks = c(0.5, 1))+format_info+facet_grid(attr~clim_class)
ggsave(filename = paste0(plot_path, "CAMELS_dist_comparison_clim1.pdf"), dpi = 'retina', width = 7, height =10)
#soils
plotdf_temp = do.call(rbind, dist_list_1py[c(14:24)])
ggplot(plotdf_temp, aes(x, colour = ggg))+stat_ecdf(size = 0.8)+col_info+xlab_info+ylab_info+theme_bw()+scale_x_continuous(breaks = c(0.5, 1))+format_info+facet_grid(attr~clim_class)
ggsave(filename = paste0(plot_path, "CAMELS_dist_comparison_soils1.pdf"), dpi = 'retina', width = 7, height =10)
#vegetation
plotdf_temp = do.call(rbind, dist_list_1py[c(25:29)])
ggplot(plotdf_temp, aes(x, colour = ggg))+stat_ecdf(size = 0.8)+col_info+xlab_info+ylab_info+theme_bw()+scale_x_continuous(breaks = c(0.5, 1))+format_info+facet_grid(attr~clim_class)
ggsave(filename = paste0(plot_path, "CAMELS_dist_comparison_vege1.pdf"), dpi = 'retina', width = 7, height =7)
#geology
plotdf_temp = do.call(rbind, dist_list_1py[c(30:32)])
ggplot(plotdf_temp, aes(x, colour = ggg))+stat_ecdf(size = 0.8)+col_info+xlab_info+ylab_info+theme_bw()+scale_x_continuous(breaks = c(0.5, 1))+format_info+facet_grid(attr~clim_class)
ggsave(filename = paste0(plot_path, "CAMELS_dist_comparison_geol1.pdf"), dpi = 'retina', width = 7, height =7)


#write plots POT3py to files ----
#loco topo
plotdf_temp = do.call(rbind, dist_list_3py[c(1:4)])
ggplot(plotdf_temp, aes(x, colour = ggg))+stat_ecdf(size = 0.8)+col_info+xlab_info+ylab_info+theme_bw()+scale_x_continuous(breaks = c(0.5, 1))+format_info+strip_theme2+facet_grid(attr~clim_class, labeller = labeller(attr = attr_fulllabel[c(1:4)]))
ggsave(filename = paste0(plot_path, "CAMELS_dist_comparison_locotopo3.pdf"), dpi = 'retina', width = 7, height =7)
#clim
plotdf_temp = do.call(rbind, dist_list_3py[c(5:13)])
ggplot(plotdf_temp, aes(x, colour = ggg))+stat_ecdf(size = 0.8)+col_info+xlab_info+ylab_info+theme_bw()+scale_x_continuous(breaks = c(0.5, 1))+format_info+strip_theme2+facet_grid(attr~clim_class, labeller = labeller(attr = attr_fulllabel[c(5:13)]))
ggsave(filename = paste0(plot_path, "CAMELS_dist_comparison_clim3.pdf"), dpi = 'retina', width = 7, height =10)
#soils
plotdf_temp = do.call(rbind, dist_list_3py[c(14:24)])
ggplot(plotdf_temp, aes(x, colour = ggg))+stat_ecdf(size = 0.8)+col_info+xlab_info+ylab_info+theme_bw()+scale_x_continuous(breaks = c(0.5, 1))+format_info+strip_theme2+facet_grid(attr~clim_class, labeller = labeller(attr = attr_fulllabel[c(14:24)]))
ggsave(filename = paste0(plot_path, "CAMELS_dist_comparison_soils3.pdf"), dpi = 'retina', width = 7, height =10)
#vegetation
plotdf_temp = do.call(rbind, dist_list_3py[c(25:29)])
ggplot(plotdf_temp, aes(x, colour = ggg))+stat_ecdf(size = 0.8)+col_info+xlab_info+ylab_info+theme_bw()+scale_x_continuous(breaks = c(0.5, 1))+format_info+strip_theme2+facet_grid(attr~clim_class, labeller = labeller(attr = attr_fulllabel[c(25:29)]))
ggsave(filename = paste0(plot_path, "CAMELS_dist_comparison_vege3.pdf"), dpi = 'retina', width = 7, height =7)
#geology
plotdf_temp = do.call(rbind, dist_list_3py[c(30:32)])
ggplot(plotdf_temp, aes(x, colour = ggg))+stat_ecdf(size = 0.8)+col_info+xlab_info+ylab_info+theme_bw()+scale_x_continuous(breaks = c(0.5, 1))+format_info+strip_theme2+facet_grid(attr~clim_class, labeller = labeller(attr = attr_fulllabel[c(30:32)]))
ggsave(filename = paste0(plot_path, "CAMELS_dist_comparison_geol3.pdf"), dpi = 'retina', width = 7, height =7)

# Plot difference in probability density distribution ---------------------------------------------------

#size dependent on event count
ggplot(diff_df_3py, aes(process, attr, col = mean_val, size = event_count))+geom_point()+theme_bw()+scale_colour_gradient2(name = "Mean diff.", limits = c(-0.25, 0.25))+scale_size(range = c(1,5), limits = range(diff_df_3py$event_count), name = "# events", breaks = c(500, 5000, 10000))+theme(axis.text.x = element_text(angle = 45, hjust = 1),axis.title = element_blank(), panel.border = element_rect(colour = 'grey80'))+strip_theme+scale_x_discrete(labels = process_fulllabel)+facet_grid(char_group~clim, scales = "free_y", space = "free_y")

ggsave(filename = paste0(plot_path, "CAMELS_meandiff_ecdf.pdf"), dpi = 'retina', width = 7, height =7)

#not size dependent
p_diff = ggplot(diff_df_3py, aes(process, attr, col = mean_val))+geom_point(size = 5)+theme_bw()+scale_colour_gradient2(name = "Mean diff.", limits = c(-0.25, 0.25), high = '#08306B')+theme(axis.text.x = element_text(angle = 45, hjust = 1),axis.title = element_blank(), panel.border = element_rect(colour = 'grey80'))+strip_theme+scale_x_discrete(labels = process_fulllabel)+facet_grid(char_group~clim, scales = "free_y", space = "free_y")
p_diff
ggsave(filename = paste0(plot_path, "CAMELS_meandiff_ecdf_equalsize.pdf"), dpi = 'retina', width = 7, height =7)

#For POT1py
ggplot(diff_df_1py, aes(process, attr, col = mean_val, size = event_count))+geom_point()+theme_bw()+scale_colour_gradient2(name = "Mean diff.", limits = c(-0.3, 0.3))+scale_size(range = c(1,5), limits = range(diff_df_1py$event_count), name = "# events", breaks = c(500, 2500, 5000))+theme(axis.text.x = element_text(angle = 45, hjust = 1),axis.title = element_blank(), panel.border = element_rect(colour = 'grey80'))+strip_theme+scale_x_discrete(labels = process_fulllabel)+facet_grid(char_group~clim, scales = "free_y", space = "free_y")

ggsave(filename = paste0(plot_path, "CAMELS_meandiff_ecdf1py.pdf"), dpi = 'retina', width = 7, height =7)


# Example distribution plot ---------------------------------------------------
#Three plots to explain steps from distribution to summarised value

#Plot only two attributes for distribution plot
plotdf_temp = do.call(rbind, dist_list_3py[c(5:13)])
example_attr = grep("p_mean|pet_mean",  plotdf_temp$attr)
plotdf_temp = plotdf_temp[example_attr,]
p1_pmean = ggplot(plotdf_temp[plotdf_temp$attr=='p_mean',], aes(x, colour = ggg))+stat_ecdf(size = 0.8)+col_info+xlab("Normalised P_mean")+ylab_info+theme_bw()+scale_x_continuous(breaks = c(0, 0.5, 1))+format_info+facet_wrap(~clim_class)+theme(legend.position = "none", panel.border = element_blank(), panel.grid.minor = element_blank())+strip_theme
p1_petmean = ggplot(plotdf_temp[plotdf_temp$attr=='pet_mean',], aes(x, colour = ggg))+stat_ecdf(size = 0.8)+col_info+xlab("Normalised PET_mean")+ylab_info+theme_bw()+scale_x_continuous(breaks = c(0, 0.5, 1))+format_info+facet_wrap(~clim_class)+theme(legend.position = "right", panel.border = element_blank(), panel.grid.minor = element_blank())+strip_theme
p1 = p1_pmean + p1_petmean+plot_layout(ncol=2, guides = "collect") & theme(legend.position = 'bottom')

plotdf_temp2 = diff_df_3py
example_attr = grep("p_mean|pet_mean",  plotdf_temp2$attr)
plotdf_temp2 = plotdf_temp2[example_attr,]

p2 = ggplot(plotdf_temp2, aes(process, attr, col = mean_val, size = event_count))+geom_point()+theme_bw()+scale_color_gradient2(name = "Mean diff.", limits = c(-0.27, 0.27), high = '#08306B')+scale_size(range = c(1,10), limits = range(diff_df_3py$event_count), name = "# events", breaks = c(500, 5000, 10000))+scale_x_discrete(labels = process_fulllabel)+theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.y = element_blank(), panel.border = element_blank())+facet_wrap(~clim)+theme(legend.position = "right",panel.spacing = unit(2, "lines"))+strip_theme+xlab("Flood Process") 

#Calculate summarised values
p_mean_soilsat = subset(plotdf_temp, attr == "p_mean" & clim_class == "Snow" & processname == "soilsat")
p_mean_longrainfall = subset(plotdf_temp, attr == "p_mean" & clim_class == "Snow" & processname == "longrainfall")
p_mean_All = subset(plotdf_temp, attr == "p_mean" & clim_class == "Snow" & processname == "All")

ecdf1 = ecdf(p_mean_soilsat$x)
ecdf2 = ecdf(p_mean_longrainfall$x)
ecdfAll = ecdf(attr_norm$p_mean[clim_index == 'Snow'])

plotdf = data.frame(X = seq(0,1, 0.01), soilsat = ecdf1(seq(0,1, 0.01)), longrainfall =  ecdf2(seq(0,1, 0.01)), All =  ecdfAll(seq(0,1, 0.01)))
plotdf_melt = data.frame(melt(plotdf, id = "X"), plotdf[rep(1:nrow(plotdf), 3),])

p3 = ggplot(plotdf_melt, aes(X, value, col = variable))+xlab("Normalised P_mean")+geom_ribbon(aes(ymin = soilsat, ymax = All), fill = muted("blue"), alpha = 0.1, show.legend = F)+geom_ribbon(aes(ymin = All, ymax = longrainfall), fill = muted("red"), alpha = 0.1, show.legend = F)+geom_line(size = 1.5)+scale_colour_manual(name = "Flood process", labels = c("Excess rainfall", "Long rainfall",  "All" ), values = c(Excessrain_col, Longrainfall_col,'black'))+annotate(geom = "text", x = 0.53, y = 0.93, label = "Negative")+annotate(geom = "text", x = 0.53, y = 0.75, label = "Positive")+theme_bw()+theme(legend.position = "none", axis.text = element_text(size = 9),  panel.border = element_blank(), panel.grid.minor = element_blank())+ylab("ECDF")
#combine plots
p1_pmean+p1_petmean+(p3|plot_spacer())+p2+plot_layout(nrow = 2)+plot_annotation(tag_levels = 'A')

ggsave(filename = paste0(plot_path, "CAMELS_example_dist.pdf"), dpi = 'retina', width = 14, height =6)

# Plot random forest predicted values ---------------------------------------------------

p_pred = ggplot(clim_pred_df, aes(100*Obs, 100*Pred, col = clim))+geom_point(shape = 1)+scale_colour_manual(limits = clim_label_vec,  values = clim_col_vec, name = '')+xlab("Observed [%]")+ylab("Prediction from cross-validation [%]")+theme_bw()+geom_abline(intercept = 0, slope =1)+theme(legend.text=element_text(size=12), panel.border = element_blank(), panel.grid.minor = element_blank())+strip_theme+facet_wrap(~process, ncol = 5,labeller = labeller(process = process_fulllabel))
p_pred
ggsave(filename = paste0(plot_path, "CAMELS_RF_predictCV.pdf"), dpi = 'retina', width = 12, height =4)


# Prediction accuracy ---------------------------------------------------

pred_Rsquared = do.call(rbind, lapply(clim_label_vec, function(clim){
  temp_out = do.call(rbind, lapply(level_vec_class, function(process){
    temp_preddf = clim_pred_df[clim_pred_df$clim == clim,]
    temp_preddf = temp_preddf[temp_preddf$process == process,]
    r2_val = r2func(temp_preddf$Obs, temp_preddf$Pred)
    return(data.frame(r2_val = r2_val, process = process, clim = clim))
  }))
  return(temp_out)
}))

#pred_Rsquared = ablation_best_full[ablation_best_full$leng==1,]

ggplot(pred_Rsquared, aes(process, r2_val))+geom_bar(stat = 'identity')+facet_wrap(~clim)
p_pred_bar = ggplot(pred_Rsquared, aes(clim, r2_val, fill = clim))+geom_bar(stat = 'identity')+scale_fill_manual(limits = clim_label_vec,  values = clim_col_vec, name = '')+ylim(0,1)+theme_bw()+theme(panel.border = element_blank(), panel.grid.minor  = element_blank(), panel.grid.major.x = element_blank(), axis.title.x = element_blank(), legend.position = 'none')+strip_theme+facet_wrap(~process,ncol = 5, labeller = labeller(process = process_fulllabel))+ylab('R-squared')

p_pred/p_pred_bar + plot_layout(guides = 'collect')+plot_annotation(tag_levels = 'A')
ggsave(filename = paste0(plot_path, "CAMELS_RF_predictCV.pdf"), dpi = 'retina', width = 12, height =6)


# Map errors in random forest predicted values ---------------------------------------------------

point_data = geom_point(data = clim_pred_df_merge, aes(x=gauge_lon, y=gauge_lat, col = pred_diff),na.rm = TRUE, size = 3, shape = 1)
p = USmap+background+point_data+scale_color_gradientn(colours=rev(pals::parula(100)), name = 'Squared prediction error',na.value = "white",limits = c(0,max(clim_pred_df_merge$pred_diff)), breaks = c(0, 0.05, 0.10, 0.15))+coord_map()+strip_theme+facet_wrap(~process, labeller = labeller(process = process_fulllabel))+guides(col = guide_colourbar(barwidth = 12.5, title.position = 'top')) #col = guide_colourbar(barwidth = 5), 
p = shift_legend2(p)

ggsave(filename = paste0(plot_path, "CAMELS_RF_squarederror_map.pdf"),plot = p, dpi = 'retina', width = 10, height =4)



# Plot summarised accumulated local effects ---------------------------------------------------

cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
names(cbPalette) <- c("Topo" , "Clim" , "Vege", "Soils" ,"Geol")
col_scale  <- scale_fill_manual(name = "remove_group",values = cbPalette)
scaleFUN <- function(x) sprintf("%.3f", x)
num_scale = scale_y_continuous(labels=scaleFUN)

#Barplot
def_theme = theme(axis.text.x = element_text(angle = 45, hjust = 1),  axis.title.x = element_blank())+strip_theme
ggplot(ALEdfagg, aes(attr, val, fill = group))+geom_col()+theme_bw()+col_scale+def_theme+num_scale+facet_grid(process~clim, scales = 'free_y', labeller = labeller(process = process_fulllabel))+ylab("Normalised Accumulated Local Effects Summary")
ggsave(filename = paste0(plot_path, "CAMELS_ALEsummary.pdf"), dpi = 'retina', width = 12, height =10)

#Point plot
def_theme = theme(axis.text.x = element_text(angle = 45, hjust = 1),  axis.title = element_blank(), panel.border = element_rect(colour = 'grey80'))+strip_theme
ALEdfagg[,"attr"] = factor(ALEdfagg$attr, levels = rev(unique(ALEdfagg$attr)))

#points scaled by no. events
ggplot(ALEdfagg, aes(process, attr, col = val, size = val))+geom_point()+scale_color_gradientn(colours=pals::brewer.blues(100), name = 'ALE')+scale_size(name = 'ALE', range = c(0, 4))+theme_bw()+def_theme+scale_x_discrete(labels = process_fulllabel)+facet_grid(group~clim, scales = 'free_y', space = 'free')
ggsave(filename = paste0(plot_path, "CAMELS_ALEsummary_points.pdf"), dpi = 'retina', width = 7, height =7)

#points not scaled
p_ALE = ggplot(ALEdfagg, aes(process, attr, col = val))+geom_point(size = 5)+scale_color_gradientn(colours=pals::brewer.blues(100), name = 'ALE')+theme_bw()+def_theme+scale_x_discrete(labels = process_fulllabel)+facet_grid(group~clim, scales = 'free_y', space = 'free')
p_ALE
ggsave(filename = paste0(plot_path, "CAMELS_ALEsummary_points_equalsize.pdf"), dpi = 'retina', width = 7, height =7)


# ALE example plot --------------------------------------------------------
char_group = c(rep("Topo", 4), rep("Clim", 9), rep("Soils", 11), rep("Vege", 5), rep("Geol", 2))
tempdf = data.frame(process = count_mechdf[,'snowmelt'], attr_cont[,-31])
tempdf = tempdf[clim_index=='Wet',]
rf = randomForest(process ~ ., data = tempdf)
X =  tempdf[which(names(tempdf) != "process")]
model = Predictor$new(rf, data = X, y = tempdf$process)

ALEdf_ex = do.call(rbind, lapply(c(1:ncol(X)), function(i){
  ALEtemp = FeatureEffect$new(model, method = 'ale', feature = colnames(X)[i])
  ALEtempdf = ALEtemp$results
  data.frame(xval = ALEtempdf[,3], fval = ALEtempdf[,1], attr = rep(colnames(X)[i], nrow(ALEtempdf)), group = rep(char_group[i], nrow(ALEtempdf)))
}))

attr_label_ex = c('p_mean'= "Mean Precipitation [mm]", 'pet_mean' = "Mean Evapotranspiration [mm]",  'water_frac' = "Water Fraction [%]")

plot_tempdf_ALE = ALEdf_ex[ALEdf_ex$attr %in% c('p_mean', 'pet_mean', 'water_frac'),]
p1 = ggplot(plot_tempdf_ALE, aes(xval, fval))+geom_bar(stat = "identity", fill = "blue", width = 0.1)+geom_line()+geom_hline(yintercept = 0, lty = 2)+theme_bw()+xlab('Attribute range')+ylab('Accumulated local effects')+theme( panel.border = element_rect(colour = 'grey80'), panel.grid = element_blank())+strip_theme+facet_wrap(~attr, scales = 'free_x', labeller = labeller(attr = attr_label_ex))

plot_tempdf = melt(tempdf[,colnames(tempdf) %in% c('p_mean', 'pet_mean', 'water_frac')])
p2 = ggplot(plot_tempdf, aes(value))+geom_histogram()+theme_bw()+xlab('Attribute range')+ylab('Frequency')+theme( panel.border = element_rect(colour = 'grey80'), panel.grid = element_blank())+strip_theme+facet_wrap(~variable, scales = 'free_x', labeller = labeller(variable = attr_label_ex))
p1+p2+plot_layout(nrow = 2)+plot_annotation(tag_levels = 'A')

ggsave(filename = paste0(plot_path, "CAMELS_example_ALE.pdf"), dpi = 'retina', width = 7, height =5)

#Example summarised values for snowmelt floods in wet catchments
ALEdfagg[ALEdfagg$attr == "p_mean",]
#0.84243064
ALEdfagg[ALEdfagg$attr == "pet_mean",]
#0.935138096
ALEdfagg[ALEdfagg$attr == "water_frac",]
#0.781624420

# Combined distribution and ALE plot ---------------------------------------------------

#requires results from ablation analysis (=R2 values for prediction with all values included)
ALE_R2 = merge(ALEdfagg, ablation_best_full[ablation_best_full$remove_group == 'All',])

ALE_R2_temp = ALE_R2
ALE_R2_temp[,'attr_label'] = rep(CAMELS_labels[-31], length.out = nrow(ALE_R2_temp))
p_ALE_r2 = ggplot(ALE_R2_temp, aes(process, attr_label, col = val, size = r2_eval_temp))+geom_point()+scale_color_gradientn(colours=pals::brewer.blues(100), name = 'ALE')+scale_size(name = 'R2',range = c(1,5))+theme_bw()+def_theme+ylab("")+scale_x_discrete(labels = process_fulllabel)+facet_grid(group~clim, scales = 'free_y', space = 'free')+ggtitle('Summarised accumulated local affects')
p_ALE_r2

diff_df_3py_temp = diff_df_3py
diff_df_3py_temp[,'attr_label'] = rep(CAMELS_labels, each = nrow(diff_df_3py_temp)/length(CAMELS_labels))

p_diff_event = ggplot(diff_df_3py_temp, aes(process, attr_label, col = mean_val, size = event_count))+geom_point()+theme_bw()+scale_colour_gradient2(name = "Mean diff.", limits = c(-0.25, 0.25))+scale_size(range = c(1,5), limits = range(diff_df_3py$event_count), name = "# events", breaks = c(500, 5000, 10000))+theme(axis.text.x = element_text(angle = 45, hjust = 1),axis.title = element_blank(), panel.border = element_rect(colour = 'grey80'))+strip_theme+scale_x_discrete(labels = process_fulllabel)+facet_grid(char_group~clim, scales = "free_y", space = "free_y")+ggtitle('Mean difference in probability distribution')
p_diff_event + p_ALE_r2 + plot_layout()+ plot_annotation(tag_levels = 'A')
ggsave(filename = paste0(plot_path, "CAMELS_importance_both_unequalsize.pdf"), dpi = 'retina', width = 12, height =7)


#POT1 combined plot -----
#requires results from ablation analysis (=R2 values for prediction with all values included)
ALEdfagg_POT1[,"attr"] = factor(ALEdfagg_POT1$attr, levels = rev(unique(ALEdfagg_POT1$attr)))

ALE_R2_POT1 = merge(ALEdfagg_POT1, ablation_best_full_POT1[ablation_best_full_POT1$remove_group == 'All',])
ALE_R2_POT1[,'attr_label'] = rep(CAMELS_labels[-31], length.out = nrow(ALE_R2_POT1))

p_ALE_r2_POT1 = ggplot(ALE_R2_POT1, aes(process, attr_label, col = val, size = r2_eval_temp))+geom_point()+scale_color_gradientn(colours=pals::brewer.blues(100), name = 'ALE')+scale_size(name = 'R2',range = c(1,5))+theme_bw()+def_theme+scale_x_discrete(labels = process_fulllabel)+ylab("")+facet_grid(group~clim, scales = 'free_y', space = 'free')+ggtitle('Summarised accumulated local affects')
p_ALE_r2_POT1


diff_df_1py_temp = diff_df_1py
diff_df_1py_temp[,'attr_label'] = rep(CAMELS_labels, each = nrow(diff_df_1py_temp)/length(CAMELS_labels))

p_diff_event_POT1 = ggplot(diff_df_1py_temp, aes(process, attr_label, col = mean_val, size = event_count))+geom_point()+theme_bw()+scale_colour_gradient2(name = "Mean diff.", limits = c(-0.25, 0.25))+scale_size(range = c(1,5), limits = range(diff_df_1py$event_count), name = "# events", breaks = c(500, 2500, 5000))+theme(axis.text.x = element_text(angle = 45, hjust = 1),axis.title = element_blank(), panel.border = element_rect(colour = 'grey80'))+strip_theme+scale_x_discrete(labels = process_fulllabel)+facet_grid(char_group~clim, scales = "free_y", space = "free_y")+ggtitle('Mean difference in probability distribution')
p_diff_event_POT1 + p_ALE_r2_POT1 + plot_layout()+ plot_annotation(tag_levels = 'A')
ggsave(filename = paste0(plot_path, "CAMELS_importance_both_unequalsize_POT1.pdf"), dpi = 'retina', width = 12, height =7)



p_diff_event + p_ALE_r2 + p_diff_event_POT1 + p_ALE_r2_POT1 + plot_layout(nrow = 2)+ plot_annotation(tag_levels = 'A')
ggsave(filename = paste0(plot_path, "CAMELS_importance_both_unequalsize_POT3_1.pdf"), dpi = 'retina', width = 12, height =12)



# Plot ablation analysis ---------------------------------------------------

cbPalette <- c('Grey30', "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
names(cbPalette) <- c('All', "Topo" , "Clim" , "Vege", "Soils" ,"Geol")
colScale <- scale_colour_manual(name = "remove_group",values = cbPalette)
p1 = ggplot(ablation_best_full, aes(6-leng, r2_eval_temp, fill = remove_group))+geom_col()+scale_x_reverse()+scale_fill_manual(values=cbPalette, name = 'Drop', drop = F)+facet_grid(process~clim,scales = 'free_x', labeller = labeller(process = process_fulllabel))+theme_bw()+theme(panel.border = element_rect(colour = 'grey80'), panel.grid = element_blank())+strip_theme+xlab('No. of attribute groups')+ylab('R2')+ggtitle('Ablation')
p2 = ggplot(ablation_rev_full, aes(6-leng, r2_eval_temp, fill = remove_group))+geom_col()+scale_x_reverse()+scale_fill_manual(values=cbPalette, name = 'Drop', drop = F)+facet_grid(process~clim,scales = 'free_x', labeller = labeller(process = process_fulllabel))+theme_bw()+theme( panel.border = element_rect(colour = 'grey80'),panel.grid = element_blank())+strip_theme+xlab('No. of attribute groups')+ylab('R2')+ggtitle('Reverse Ablation')

p1 + p2 +  plot_layout(guides = 'collect')+ plot_annotation(tag_levels = 'A')
ggsave(filename = paste0(plot_path, "CAMELS_Ablation.pdf"), dpi = 'retina', width = 12, height =8)







