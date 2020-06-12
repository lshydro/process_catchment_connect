###########################################################-
# Objective: Function file to calculate flood generating
#mechanisms
# Author: Lina Stein, University of Bristol
# note: 
###########################################################-

#variable_method_flood_mech_dataseries.R

#
#  max. daily rainfall  -----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Calc_Pdaily_dt = function(MSWEP_cat_vec){
  temp_df = data.frame(ind_years, doY, MSWEP_cat_vec)
  colnames(temp_df) = c("Year", "doY", "MaxPdaily")
  dt <- data.table(temp_df)
  dtm = dt[, .SD[which.max(MaxPdaily),], by=Year]
}

#
#  max. weekly rainfall  -----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Calc_Pweekly = function(MSWEP_cat_vec){
  av_7_day = rollapply(MSWEP_cat_vec, 7, sum, partial = T, align = "right") #partial: apply for start values for smaller range
  temp_df = data.frame(ind_years, doY, av_7_day)
  colnames(temp_df) = c("Year", "doY", "MaxPweek")
  dt <- data.table(temp_df)
  dtm = dt[, .SD[which.max(MaxPweek),], by=Year]
}


#
#  soil routine function  -----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


eff_func = function(AWC_in, ET_scale, GLEAM_mat, path_out){
  
  ind_years = years(ind_df[,1])
  doY = format(ind_df[,1], format = "%j")
  
  #effective Precipitation
  P_eff = lapply(c(1:length(MSWEP_reduced[1,])), FUN = function(x){ 
    P_vec = MSWEP_reduced[,x]
    ET_vec = GLEAM_mat[,x]
    Ssmax = AWC_in[x]
    
    #Soil storage
    Ssoil_vec = rep(NA, length(P_vec))
    Ssoil_vec[1] = 0
    #Precipitation excess
    P_exc_vec = rep(NA, length(P_vec))
    
    #only do calculations if soil value is available, otherwise P_exc_vec stays NA
    if(!is.na(Ssmax)){
      P_exc_vec[1] = 0
      for (i in c(2:length(P_vec))){
        #soil storage plus precipitation
        Ssoil_vec[i] = Ssoil_vec[i-1]+P_vec[i]
        if (Ssoil_vec[i]>Ssmax){ 
          #if soil storage exceeds maximum soil storage, excess becomes runoff
          P_exc_vec[i] = Ssoil_vec[i]-Ssmax
          #and soil storage is limited at maximum soil storage
          Ssoil_vec[i] = Ssmax
        } else {
          #if maximum soil storage is not exceeded, no runoff occurs
          P_exc_vec[i] = 0
        }
        #soil storage is reduced by (delayed) evaporation
        Ssoil_vec[i] = Ssoil_vec[i] - min(c((ET_scale * ET_vec[i]), Ssoil_vec[i]))
      }
    }
    
    P_exc_vec
  })
  
  #result is list of one long vector, turn back into matrix
  P_effmat = matrix(unlist(P_eff), nrow = nrow(MSWEP_reduced), ncol = ncol(MSWEP_reduced), byrow = F)
  save(P_effmat, file = paste0(path_out, "Peffmat.Rdata"))
  
  
  
  Calc_Peff = function(Peff_vec){
    temp_df = data.frame(ind_years, doY, Peff_vec)
    colnames(temp_df) = c("Year", "doY", "MaxPeff")
    #only calculate max value if values are available, otherwise return data table filled with NA
    if(!all(is.na(Peff_vec))){
      dt <- data.table(temp_df)
      dtm = dt[, .SD[which.max(MaxPeff),], by=Year]
    } else {
      na_df = data.frame(unique(ind_years), rep(NA, length(unique(ind_years))), rep(NA, length(unique(ind_years))))
      colnames(na_df) = c("Year", "doY", "MaxPeff")
      dtm = data.table(na_df)
    }
  }
  
  P_effmax = apply(P_effmat, 2, FUN = Calc_Peff)
  save(P_effmax, file = paste0(path_out, "Peffmax.Rdata"))
  
  # Calculate soil saturation -----
  
  soil_bucket_vec = lapply(c(1:length(MSWEP_reduced[1,])), FUN = function(x){ 
    P_vec = MSWEP_reduced[,x]
    ET_vec = GLEAM_mat[,x]
    Ssmax = AWC_in[x]
    
    #Soil storage
    Ssoil_vec = rep(NA, length(P_vec))
    Ssoil_vec[1] = 0
    #Precipitation excess
    P_exc_vec = rep(NA, length(P_vec))
    
    #only do calculations if soil value is available, otherwise P_exc_vec stays NA
    if(!is.na(Ssmax)){
      P_exc_vec[1] = 0
      for (i in c(2:length(P_vec))){
        #soil storage plus precipitation
        Ssoil_vec[i] = Ssoil_vec[i-1]+P_vec[i]
        if (Ssoil_vec[i]>Ssmax){ 
          #if soil storage exceeds maximum soil storage, excess becomes runoff
          P_exc_vec[i] = Ssoil_vec[i]-Ssmax
          #and soil storage is limited at maximum soil storage
          Ssoil_vec[i] = Ssmax
        } else {
          #if maximum soil storage is not exceeded, no runoff occurs
          P_exc_vec[i] = 0
        }
        #soil storage is reduced by (delayed) evaporation
        Ssoil_vec[i] = Ssoil_vec[i] - min(c((ET_scale * ET_vec[i]), Ssoil_vec[i]))
      }
    }
    
    Ssoil_vec
  })
  
  #result is list of one long vector, turn back into matrix
  soil_bucket = matrix(unlist(soil_bucket_vec), nrow = nrow(MSWEP_reduced), ncol = ncol(MSWEP_reduced), byrow = F)
  save(soil_bucket, file = paste0(path_out, "soil_bucket.Rdata"))
  
  soil_sat = apply(soil_bucket, 1, FUN = function(x){x/AWC_in})
  save(soil_sat, file = paste0(path_out, "soil_sat.Rdata"))

}



#
#  snowmelt routine function -----
#

snow_func = function(version_num, fdd_in, Tcrit_in){
  version = version_num
  file_output_path = paste0("C:/ls16959/Data/04_Calculations/mech_compare/", version, "_")
  fdd = fdd_in #(mm/d/k)
  ind_years = years(ind_df[,1])
  doY = format(ind_df[,1], format = "%j")
  
  P_snowmat = matrix(NA, nrow = nrow(MSWEP_reduced), ncol = ncol(MSWEP_reduced))
  
  for (j in c(1:length(MSWEP_reduced[1,]))){ 
    #print(j)
    Ta_vec = BEST_reduced[,j]
    P_vec = MSWEP_reduced[,j]
    #Snow storage
    Ssnow_vec = rep(NA, length(P_vec))
    Ssnow_vec[1] = 0
    #Snow melt
    P_snow_vec = rep(NA, length(P_vec))
    P_snow_vec[1] = 0
    
    for (i in c(2:length(P_vec))){
      #if temperature data is not available snowmelt is NA
      if (is.na(Ta_vec[i])){
        P_snow_vec[i] = NA
      } 
      else {
        #if critical temperature data is not available snowmelt is NA
        if(is.na(Tcrit_in[j])){
          P_snow_vec[i] = NA
        } else {
          if (Ta_vec[i]<Tcrit_in[j]){
            Ssnow_vec[i] = Ssnow_vec[i-1] + P_vec[i]
            P_snow_vec[i] = 0
          } 
          else{
            #if temperature above Tcrit, snowmelt happens
            #snowmelt is either calculated from degree day with the max amount of snow to melt is the snow storage
            P_snow_vec[i] = min(c(fdd*max(c(Ta_vec[i]-Tcrit_in[j], 0)), Ssnow_vec[i-1])) + P_vec[i]
            Ssnow_vec[i] = max(Ssnow_vec[i-1] - P_snow_vec[i], 0)
            #if there is no snow to melt, no snowmelt happens, only rainfall
            if (Ssnow_vec[i]<=0){
              P_snow_vec[i] = 0
            }
          }
        }
      }
    }
    P_snowmat[,j] = P_snow_vec
  }
  save(P_snowmat, file = paste0(file_output_path, "Psnowmat.Rdata"))
  
  #which.max, if only NA, value will be discarded
  Calc_Psnow = function(Psnow_vec){
    temp_df = data.frame(ind_years, doY, Psnow_vec)
    colnames(temp_df) = c("Year", "doY", "MaxPsnow")
    dt <- data.table(temp_df)
    dtm = dt[, .SD[which.max(MaxPsnow),], by=Year]
  }
  
  P_snowmax = apply(P_snowmat, 2, FUN = Calc_Psnow)
  save(P_snowmax, file = paste0(file_output_path, "Psnowmax.Rdata"))
  
} 

#
#  coupled soil-snow routine -----
#
#####
 
# soil_snow_func = function(AWC_in, ET_scale, ET_mat, rain_mat, Temp_mat, path_out, fdd_in, Tcrit_in, RoS = T){
#   #AWC_in: soil storage value
#   #ET_scale: 1 for actual ET, 0.75 for potential ET
#   #ET_mat: evapotranspiration matrix (col: catchments, row: days)
#   #rain_mat: precipitation matrix (col: catchments, row: days)
#   #Temp_mat: temperature matrix (col: catchments, row: days)
#   #path_out: file path to store output files
#   #fdd_in: melt rate
#   #Tcrit_in: critical temperature below which rain turns to snow
#   #RoS: output of snowroutine with P or without P
#   
#   
#   fdd = fdd_in #(mm/d/k)
#   
#   P_snowmat = matrix(NA, nrow = nrow(rain_mat), ncol = ncol(rain_mat))
#   P_effmat = matrix(NA, nrow = nrow(rain_mat), ncol = ncol(rain_mat))
#   P_satmat = matrix(NA, nrow = nrow(rain_mat), ncol = ncol(rain_mat))
#   P_snowstrgmat = matrix(NA, nrow = nrow(rain_mat), ncol = ncol(rain_mat))
#   
#   for (j in c(1:length(rain_mat[1,]))){ 
#     Ta_vec = Temp_mat[,j]
#     P_vec = rain_mat[,j]
#     ET_vec = ET_mat[,j] * ET_scale
#     ET_vec[ET_vec<0] = 0
#     Ssmax = AWC_in[j]
#     Tcrit = Tcrit_in[j]
#     
#     #Snow storage
#     Ssnow_vec = rep(NA, length(P_vec))
#     Ssnow_vec[1] = 0
#     #Snow melt
#     P_snow_vec = rep(NA, length(P_vec))
#     P_snow_vec[1] = 0
#     #Soil storage
#     Ssoil_vec = rep(NA, length(P_vec))
#     Ssoil_vec[1] = 0
#     #Soil overflow
#     soil_overflw = rep(NA, length(P_vec))
#     
#     for (i in c(2:length(P_vec))){
#       # ~~~~~ calculate snow storage ~~~~~ 
#       if (is.na(Tcrit)){
#         #if critical temperature data is not available snowstorage is NA
#         Ssnow_vec[i] = NA
#       } else if(is.na(Ta_vec[i])){
#         #if temperature data is not available snowstorage is NA
#         Ssnow_vec[i] = NA
#       } else if(Ta_vec[i]<Tcrit){
#         Ssnow_vec[i] = Ssnow_vec[i-1] + P_vec[i]
#       } else {
#         snowmelt_timestep = min(c(fdd*max(c(Ta_vec[i]-Tcrit_in[j], 0)), Ssnow_vec[i-1]))
#         Ssnow_vec[i] = Ssnow_vec[i-1]-snowmelt_timestep
#       }
#       
#       # ~~~~~ calculate soil storage ~~~~~
#       #relevant for next step
# 
#       if (is.na(Tcrit)){
#         #if critical temperature data is not available soil stroage is filling independent of temperature
#         Ssoil_vec[i] = ifelse(is.na(Ssmax), NA, max(min(Ssoil_vec[i-1]+P_vec[i], Ssmax) - ET_vec[i],0)) #range of soil storage: 0, Ssmax
#       } else if(is.na(Ta_vec[i])){
#         #if temperature data is not available soil storage is filling independent of temperature
#         Ssoil_vec[i] = ifelse(is.na(Ssmax), NA, max(min(Ssoil_vec[i-1]+P_vec[i], Ssmax) - ET_vec[i],0)) #range of soil storage: 0, Ssmax
#       } else if(Ta_vec[i]<Tcrit){
#         #if T below Tcrit soil filling/ET is paused
#         Ssoil_vec[i] = Ssoil_vec[i-1]
#       } else if(Ssnow_vec[i] == 0 | is.na(Ssnow_vec[i])){
#         #soil filling only when there is no snow cover or snow cover is NA
#         Ssoil_vec[i] = ifelse(is.na(Ssmax), NA, max(min(Ssoil_vec[i-1]+P_vec[i], Ssmax) - ET_vec[i],0)) 
#       } else {
#         #during snow cover, soil filling/ET is paused
#         Ssoil_vec[i] = Ssoil_vec[i-1]
#       }
#       
#       # ~~~~~ calculate soil storage overflow ~~~~~
#       if (is.na(Tcrit)){
#         #if critical temperature data is not available soil overflow is happening independent of temperature
#         soil_overflw[i] = ifelse(is.na(Ssmax), NA, max(Ssoil_vec[i-1]+P_vec[i] - Ssmax, 0))
#       } else if(is.na(Ta_vec[i])){
#         #if temperature data is not available soil overflow is happening independent of temperature
#         soil_overflw[i] = ifelse(is.na(Ssmax), NA, max(Ssoil_vec[i-1]+P_vec[i] - Ssmax, 0))
#       } else if(Ta_vec[i]<Tcrit){
#         soil_overflw[i] = 0
#       } else if(Ssnow_vec[i] == 0 | is.na(Ssnow_vec[i])){
#         #soil overflow only when there is no snow cover or snow cover is NA
#         soil_overflw[i] = ifelse(is.na(Ssmax), NA, max(Ssoil_vec[i-1]+P_vec[i] - Ssmax, 0))
#       } else {
#         #during snow cover, soil overflow is zero
#         soil_overflw[i] = 0
#       }
#       
#       # ~~~~~ calculate snowmelt ~~~~~
#       if (is.na(Tcrit)){
#         #if critical temperature data is not available snowstorage is NA
#         P_snow_vec[i] = NA
#       } else if(is.na(Ta_vec[i])){
#         #if temperature data is not available snowstorage is NA
#         P_snow_vec[i] = NA
#       } else if(Ta_vec[i]<Tcrit){
#         P_snow_vec[i] = 0
#       } else if(Ssnow_vec[i-1]>0) {
#         if(RoS){
#           P_snow_vec[i] = min(c(fdd*max(c(Ta_vec[i]-Tcrit_in[j], 0)), Ssnow_vec[i-1]))+P_vec[i]
#         } else {
#           P_snow_vec[i] = min(c(fdd*max(c(Ta_vec[i]-Tcrit_in[j], 0)), Ssnow_vec[i-1]))
#         }
#         
#       } else {
#         P_snow_vec[i] = 0
#       }
#     }
#     P_effmat[,j] = soil_overflw
#     P_satmat[,j] = Ssoil_vec/Ssmax
#     P_snowmat[,j] = P_snow_vec
#     P_snowstrgmat[,j] = Ssnow_vec
#     
#     
#   }
#   save(P_effmat, file = paste0(path_out, "Peffmat.Rdata"))
#   save(P_satmat, file = paste0(path_out, "Psatmat.Rdata"))
#   save(P_snowmat, file = paste0(path_out, "Psnowmat.Rdata"))
#   save(P_snowstrgmat, file = paste0(path_out, "Psnowstrgmat.Rdata"))
#   
#   
# } 
# 
# 




# soil_snow_func = function(AWC_in, ET_scale, ET_mat, rain_mat, Temp_mat, path_out, fdd_in, Tcrit_in, RoS = T, date_vec){
#   #AWC_in: soil storage value
#   #ET_scale: 1 for actual ET, 0.75 for potential ET
#   #ET_mat: evapotranspiration matrix (col: catchments, row: days)
#   #rain_mat: precipitation matrix (col: catchments, row: days)
#   #Temp_mat: temperature matrix (col: catchments, row: days)
#   #path_out: file path to store output files
#   #fdd_in: melt rate
#   #Tcrit_in: critical temperature below which rain turns to snow
#   #RoS: output of snowroutine with P or without P
#   
#   
#   fdd = fdd_in #(mm/d/k)
#   
#   P_snowmat = matrix(NA, nrow = nrow(rain_mat), ncol = ncol(rain_mat))
#   P_effmat = matrix(NA, nrow = nrow(rain_mat), ncol = ncol(rain_mat))
#   P_satmat = matrix(NA, nrow = nrow(rain_mat), ncol = ncol(rain_mat))
#   P_snowstrgmat = matrix(NA, nrow = nrow(rain_mat), ncol = ncol(rain_mat))
#   
#   for (j in c(1:length(rain_mat[1,]))){ 
#     Ta_vec = Temp_mat[,j]
#     P_vec = rain_mat[,j]
#     ET_vec = ET_mat[,j] * ET_scale
#     ET_vec[ET_vec<0] = 0
#     Ssmax = AWC_in[j]
#     Tcrit = Tcrit_in[j]
#     #set snow storage to zero on specific date in summer (Freudiger et al. 2014)
#     
#     #option 1: end of hottest month on average (use 28th of month for potential February dates)
#     month_ind = format(date_vec, "%m")
#     av_month_temp = aggregate(Ta_vec, by = list(month_ind), mean, na.rm = T)[,2]
#     max_temp_month = which.max(av_month_temp)
#     zero_date = format(as.Date(paste0(max_temp_month, "-28"), "%m-%d"), "%m-%d") #convert to date but then cut year
#     zero_ind = format(date_vec,"%m-%d") == zero_date
#     
#     # #option 2: fixed date by hemisphere
#     # #1 for north (2nd of August), 2 for south (2nd of February) 
#     # if(hemisph_ind[j] == 1){
#     #   zero_date = format(as.Date("08-02", "%m-%d"), "%m-%d")
#     #   zero_ind = format(date_vec,"%m-%d") == zero_date
#     # } else if(hemisph_ind[j] == 2){
#     #   zero_date = format(as.Date("02-02", "%m-%d"), "%m-%d")
#     #   zero_ind = format(date_vec,"%m-%d") == zero_date
#     # }
#     
#     #Snow storage
#     Ssnow_vec = rep(NA, length(P_vec))
#     Ssnow_vec[1] = 0
#     #Snow melt
#     P_snow_vec = rep(NA, length(P_vec))
#     P_snow_vec[1] = 0
#     #Soil storage
#     Ssoil_vec = rep(NA, length(P_vec))
#     Ssoil_vec[1] = 0
#     #Soil overflow
#     soil_overflw = rep(NA, length(P_vec))
#     
#     for (i in c(2:length(P_vec))){
#       # ~~~~~ calculate snow storage ~~~~~ 
#       if (is.na(Tcrit)){
#         #if critical temperature data is not available snowstorage is NA
#         Ssnow_vec[i] = NA
#       } else if(is.na(Ta_vec[i])){
#         #if temperature data is not available snowstorage is NA
#         Ssnow_vec[i] = NA
#       } else if(Ta_vec[i]<Tcrit){
#         Ssnow_vec[i] = Ssnow_vec[i-1] + P_vec[i]
#         if(zero_ind[i]){
#           Ssnow_vec[i] = 0
#         }
#       } else {
#         snowmelt_timestep = min(c(fdd*max(c(Ta_vec[i]-Tcrit_in[j], 0)), Ssnow_vec[i-1]))
#         Ssnow_vec[i] = Ssnow_vec[i-1]-snowmelt_timestep
#         if(zero_ind[i]){
#           Ssnow_vec[i] = 0
#         }
#       }
#       
#       # ~~~~~ calculate soil storage ~~~~~
#       #relevant for next step
#       
#       if (is.na(Tcrit)){
#         #if critical temperature data is not available soil stroage is filling independent of temperature
#         Ssoil_vec[i] = ifelse(is.na(Ssmax), NA, max(min(Ssoil_vec[i-1]+P_vec[i], Ssmax) - ET_vec[i],0)) #range of soil storage: 0, Ssmax
#       } else if(is.na(Ta_vec[i])){
#         #if temperature data is not available soil storage is filling independent of temperature
#         Ssoil_vec[i] = ifelse(is.na(Ssmax), NA, max(min(Ssoil_vec[i-1]+P_vec[i], Ssmax) - ET_vec[i],0)) #range of soil storage: 0, Ssmax
#       } else if(Ta_vec[i]<Tcrit){
#         #if T below Tcrit soil filling/ET is paused
#         Ssoil_vec[i] = Ssoil_vec[i-1]
#       } else if(Ssnow_vec[i] == 0 | is.na(Ssnow_vec[i])){
#         #soil filling only when there is no snow cover or snow cover is NA
#         Ssoil_vec[i] = ifelse(is.na(Ssmax), NA, max(min(Ssoil_vec[i-1]+P_vec[i], Ssmax) - ET_vec[i],0)) 
#       } else {
#         #during snow cover, soil filling/ET is paused
#         Ssoil_vec[i] = Ssoil_vec[i-1]
#       }
#       
#       # ~~~~~ calculate soil storage overflow ~~~~~
#       if (is.na(Tcrit)){
#         #if critical temperature data is not available soil overflow is happening independent of temperature
#         soil_overflw[i] = ifelse(is.na(Ssmax), NA, max(Ssoil_vec[i-1]+P_vec[i] - Ssmax, 0))
#       } else if(is.na(Ta_vec[i])){
#         #if temperature data is not available soil overflow is happening independent of temperature
#         soil_overflw[i] = ifelse(is.na(Ssmax), NA, max(Ssoil_vec[i-1]+P_vec[i] - Ssmax, 0))
#       } else if(Ta_vec[i]<Tcrit){
#         soil_overflw[i] = 0
#       } else if(Ssnow_vec[i] == 0 | is.na(Ssnow_vec[i])){
#         #soil overflow only when there is no snow cover or snow cover is NA
#         soil_overflw[i] = ifelse(is.na(Ssmax), NA, max(Ssoil_vec[i-1]+P_vec[i] - Ssmax, 0))
#       } else {
#         #during snow cover, soil overflow is zero
#         soil_overflw[i] = 0
#       }
#       
#       # ~~~~~ calculate snowmelt ~~~~~
#       if (is.na(Tcrit)){
#         #if critical temperature data is not available snowstorage is NA
#         P_snow_vec[i] = NA
#       } else if(is.na(Ta_vec[i])){
#         #if temperature data is not available snowstorage is NA
#         P_snow_vec[i] = NA
#       } else if(Ta_vec[i]<Tcrit){
#         P_snow_vec[i] = 0
#       } else if(Ssnow_vec[i-1]>0) {
#         if(RoS){
#           P_snow_vec[i] = min(c(fdd*max(c(Ta_vec[i]-Tcrit_in[j], 0)), Ssnow_vec[i-1]))+P_vec[i]
#         } else {
#           P_snow_vec[i] = min(c(fdd*max(c(Ta_vec[i]-Tcrit_in[j], 0)), Ssnow_vec[i-1]))
#         }
#         
#       } else {
#         P_snow_vec[i] = 0
#       }
#     }
#     P_effmat[,j] = soil_overflw
#     P_satmat[,j] = Ssoil_vec/Ssmax
#     P_snowmat[,j] = P_snow_vec
#     P_snowstrgmat[,j] = Ssnow_vec
#     
#     
#   }
#   save(P_effmat, file = paste0(path_out, "Peffmat.Rdata"))
#   save(P_satmat, file = paste0(path_out, "Psatmat.Rdata"))
#   save(P_snowmat, file = paste0(path_out, "Psnowmat.Rdata"))
#   save(P_snowstrgmat, file = paste0(path_out, "Psnowstrgmat.Rdata"))
#   
#   
# } 



# 
# soil_snow_func = function(AWC_in, ET_scale, ET_mat, rain_mat, Temp_mat, path_out, fdd_in, Tcrit_in, RoS = T, date_vec){
#   #AWC_in: soil storage value
#   #ET_scale: 1 for actual ET, 0.75 for potential ET
#   #ET_mat: evapotranspiration matrix (col: catchments, row: days)
#   #rain_mat: precipitation matrix (col: catchments, row: days)
#   #Temp_mat: temperature matrix (col: catchments, row: days)
#   #path_out: file path to store output files
#   #fdd_in: melt rate
#   #Tcrit_in: critical temperature below which rain turns to snow
#   #RoS: output of snowroutine with P or without P
#   
#   
#   fdd = fdd_in #(mm/d/k)
#   
#   P_snowmat = matrix(NA, nrow = nrow(rain_mat), ncol = ncol(rain_mat))
#   P_effmat = matrix(NA, nrow = nrow(rain_mat), ncol = ncol(rain_mat))
#   P_satmat = matrix(NA, nrow = nrow(rain_mat), ncol = ncol(rain_mat))
#   P_snowstrgmat = matrix(NA, nrow = nrow(rain_mat), ncol = ncol(rain_mat))
#   
#   for (j in c(1:length(rain_mat[1,]))){ 
#     Ta_vec = Temp_mat[,j]
#     P_vec = rain_mat[,j]
#     ET_vec = ET_mat[,j] * ET_scale
#     ET_vec[ET_vec<0] = 0
#     Ssmax = AWC_in[j]
#     Tcrit = Tcrit_in[j]
#     #set snow storage to zero on specific date in summer (Freudiger et al. 2014)
#     
#     #option 1: end of hottest month on average (use 28th of month for potential February dates)
#     month_ind = format(date_vec, "%m")
#     av_month_temp = aggregate(Ta_vec, by = list(month_ind), mean, na.rm = T)[,2]
#     max_temp_month = which.max(av_month_temp)
#     zero_date = format(as.Date(paste0(max_temp_month, "-28"), "%m-%d"), "%m-%d") #convert to date but then cut year
#     zero_ind = format(date_vec,"%m-%d") == zero_date
#     
#     # #option 2: fixed date by hemisphere
#     # #1 for north (2nd of August), 2 for south (2nd of February) 
#     # if(hemisph_ind[j] == 1){
#     #   zero_date = format(as.Date("08-02", "%m-%d"), "%m-%d")
#     #   zero_ind = format(date_vec,"%m-%d") == zero_date
#     # } else if(hemisph_ind[j] == 2){
#     #   zero_date = format(as.Date("02-02", "%m-%d"), "%m-%d")
#     #   zero_ind = format(date_vec,"%m-%d") == zero_date
#     # }
#     
#     #Snow storage
#     Ssnow_vec = rep(NA, length(P_vec))
#     Ssnow_vec[1] = 0
#     #Snow melt
#     P_snow_vec = rep(NA, length(P_vec))
#     P_snow_vec[1] = 0
#     #Soil storage
#     Ssoil_vec = rep(NA, length(P_vec))
#     Ssoil_vec[1] = 0
#     #Soil overflow
#     soil_overflw = rep(NA, length(P_vec))
#     
#     for (i in c(2:length(P_vec))){
#       # ~~~~~ calculate snow storage ~~~~~ 
#       if (is.na(Tcrit)){
#         #if critical temperature data is not available snowstorage is NA
#         Ssnow_vec[i] = NA
#       } else if(is.na(Ta_vec[i])){
#         #if temperature data is not available snowstorage is NA
#         Ssnow_vec[i] = NA
#       } else if(Ta_vec[i]<Tcrit){
#         Ssnow_vec[i] = Ssnow_vec[i-1] + P_vec[i]
#         if(zero_ind[i]){
#           Ssnow_vec[i] = 0
#         }
#       } else {
#         snowmelt_timestep = min(c(fdd*max(c(Ta_vec[i]-Tcrit_in[j], 0)), Ssnow_vec[i-1]))
#         Ssnow_vec[i] = Ssnow_vec[i-1]-snowmelt_timestep
#         if(zero_ind[i]){
#           Ssnow_vec[i] = 0
#         }
#       }
#       
#       # ~~~~~ calculate soil storage ~~~~~
#       #relevant for next step
#       
#       if (is.na(Tcrit)){
#         #if critical temperature data is not available soil stroage is filling independent of temperature
#         Ssoil_vec[i] = ifelse(is.na(Ssmax), NA, max(min(Ssoil_vec[i-1]+P_vec[i], Ssmax) - ET_vec[i],0)) #range of soil storage: 0, Ssmax
#       } else if(is.na(Ta_vec[i])){
#         #if temperature data is not available soil storage is filling independent of temperature
#         Ssoil_vec[i] = ifelse(is.na(Ssmax), NA, max(min(Ssoil_vec[i-1]+P_vec[i], Ssmax) - ET_vec[i],0)) #range of soil storage: 0, Ssmax
#       } else if(Ta_vec[i]<Tcrit){
#         #if T below Tcrit soil filling/ET is paused
#         Ssoil_vec[i] = Ssoil_vec[i-1]
#       } else if(Ssnow_vec[i] == 0 | is.na(Ssnow_vec[i])){
#         #soil filling only when there is no snow cover or snow cover is NA
#         Ssoil_vec[i] = ifelse(is.na(Ssmax), NA, max(min(Ssoil_vec[i-1]+P_vec[i], Ssmax) - ET_vec[i],0)) 
#       } else {
#         #during snow cover, soil filling/ET is paused
#         Ssoil_vec[i] = Ssoil_vec[i-1]
#       }
#       
#       # ~~~~~ calculate soil storage overflow ~~~~~
#       if (is.na(Tcrit)){
#         #if critical temperature data is not available soil overflow is happening independent of temperature
#         soil_overflw[i] = ifelse(is.na(Ssmax), NA, max(Ssoil_vec[i-1]+P_vec[i] - Ssmax, 0))
#       } else if(is.na(Ta_vec[i])){
#         #if temperature data is not available soil overflow is happening independent of temperature
#         soil_overflw[i] = ifelse(is.na(Ssmax), NA, max(Ssoil_vec[i-1]+P_vec[i] - Ssmax, 0))
#       } else if(Ta_vec[i]<Tcrit){
#         soil_overflw[i] = 0
#       } else if(Ssnow_vec[i] == 0 | is.na(Ssnow_vec[i])){
#         #soil overflow only when there is no snow cover or snow cover is NA
#         soil_overflw[i] = ifelse(is.na(Ssmax), NA, max(Ssoil_vec[i-1]+P_vec[i] - Ssmax, 0))
#       } else {
#         #during snow cover, soil overflow is zero
#         soil_overflw[i] = 0
#       }
#       
#       # ~~~~~ calculate snowmelt ~~~~~
#       if (is.na(Tcrit)){
#         #if critical temperature data is not available snowstorage is NA
#         P_snow_vec[i] = NA
#       } else if(is.na(Ta_vec[i])){
#         #if temperature data is not available snowstorage is NA
#         P_snow_vec[i] = NA
#       } else if(Ta_vec[i]<Tcrit){
#         P_snow_vec[i] = 0
#       } else if(Ssnow_vec[i-1]>0) {
#         if(RoS){
#           P_snow_vec[i] = min(c(fdd*max(c(Ta_vec[i]-Tcrit_in[j], 0)), Ssnow_vec[i-1]))+P_vec[i]
#         } else {
#           P_snow_vec[i] = min(c(fdd*max(c(Ta_vec[i]-Tcrit_in[j], 0)), Ssnow_vec[i-1]))
#         }
#         
#       } else {
#         P_snow_vec[i] = 0
#       }
#     }
#     P_effmat[,j] = soil_overflw
#     P_satmat[,j] = Ssoil_vec/Ssmax
#     P_snowmat[,j] = P_snow_vec
#     P_snowstrgmat[,j] = Ssnow_vec
#     
#     
#   }
#   save(P_effmat, file = paste0(path_out, "Peffmat.Rdata"))
#   save(P_satmat, file = paste0(path_out, "Psatmat.Rdata"))
#   save(P_snowmat, file = paste0(path_out, "Psnowmat.Rdata"))
#   save(P_snowstrgmat, file = paste0(path_out, "Psnowstrgmat.Rdata"))
#   
#   
# } 
#####


# Currently used Soil-Snow routine ----------------------------------------


snowmelt_func = function(fdd_func, T_func, Tcrit_func, Ssnow_func){
  #snowmelt is dependent on fdd and temperature difference between T and Tcrit or snow storage, whichever is smaller
  min(fdd_func*max(T_func-Tcrit_func, 0), Ssnow_func)
}
soilstrg_func = function(Ssoil_prev_func, P_func, Ssmax_func, ET_func){
  #soil storage is previous soil storage + rainfall (limited by maximum soil storage) minus ET. Value cannot go below zero
  max(min(Ssoil_prev_func+P_func, Ssmax_func) - ET_func,0)
}
soiloverflow_func = function(Ssoil_prev_func, P_func, Ssmax_func){
  #soil overflow is previous soil storage + rainfall minus max. soil storage. Value cannot go below zero
  max(Ssoil_prev_func+P_func - Ssmax_func, 0)
}

soil_snow_func = function(AWC_in, ET_mat, rain_mat, Temp_mat, fdd_in, Tcrit_in, date_vec){
  #AWC_in: soil storage value
  #ET_mat: evapotranspiration matrix (col: catchments, row: days)
  #rain_mat: precipitation matrix (col: catchments, row: days)
  #Temp_mat: temperature matrix (col: catchments, row: days)
  #path_out: file path to store output files
  #fdd_in: melt rate
  #Tcrit_in: critical temperature below which rain turns to snow
  #date_vec: time stamp for timeseries
  
  
  fdd = fdd_in #(mm/d/k)
  month_ind = format(date_vec, "%m")
  month_day_ind = format(date_vec,"%m-%d")
  
  #if Tcrit is not available (for example in the southern hemisphere), set value to 1 degree which is the average Jennings et al found in the northern hemisphere

  
  class_input_df_list = mapply(FUN = function(P_vec, Ta_vec, ET_vec, Ssmax, Tcrit){
    ET_vec[ET_vec<0] = 0
    
    if(is.na(Tcrit)){
      Tcrit_used = 1
    } else{
      Tcrit_used = Tcrit
    }

    #Snow storage accumulation
    #end of hottest month on average (use 28th of month for potential February dates)
    
    av_month_temp = aggregate(Ta_vec, by = list(month_ind), mean, na.rm = T)[,2]
    max_temp_month = as.character(which.max(av_month_temp))
    if(nchar(max_temp_month)==1){
      zero_date = paste0("0", max_temp_month, "-28")
    } else {
      zero_date = paste0(max_temp_month, "-28")
    }
    zero_ind = month_day_ind == zero_date
    
    #Vec pre-allocation
    #Snow storage
    Ssnow_vec = rep(NA, length(P_vec))
    Ssnow_vec[1] = 0
    #Snow melt
    P_snow_vec = rep(NA, length(P_vec))
    P_snow_vec[1] = 0
    #Soil storage
    Ssoil_vec = rep(NA, length(P_vec))
    Ssoil_vec[1] = 0
    #Soil overflow
    soil_overflw = rep(NA, length(P_vec))
    
    #only do calculations if soil storage value is available
    if(!is.na(Ssmax)){
      #only caculate soil calculations if Tcrit_used = NA
      if(is.na(Tcrit_used)){
        for(i in c(2:length(Ssoil_vec))){
          Ssoil_vec[i] = soilstrg_func(Ssoil_vec[i-1], P_vec[i], Ssmax, ET_vec[i])
          soil_overflw[i] = soiloverflow_func(Ssoil_vec[i-1], P_vec[i], Ssmax)
        }
      } else {
        #calculate soil and snow routine
        for(i in c(2:length(Ssoil_vec))){
          #calculate snow routine
          #snow storage
          if(is.na(Ta_vec[i])){
            #if temperature data is not available snowstorage is NA
            Ssnow_vec[i] = NA
          } else if(Ta_vec[i]<Tcrit_used){
            Ssnow_vec[i] = Ssnow_vec[i-1] + P_vec[i]
            P_snow_vec[i] = 0
            if(zero_ind[i]){
              Ssnow_vec[i] = 0
            }
          } else {
            P_snow_vec[i] = snowmelt_func(fdd, Ta_vec[i], Tcrit_used, Ssnow_vec[i-1])
            Ssnow_vec[i] = max(Ssnow_vec[i-1]-P_snow_vec[i],0)
            if(zero_ind[i]){
              Ssnow_vec[i] = 0
            }
          }
        }
        
        #"Pause" soil routine during snow cover
        P_vec_soil = P_vec
        P_vec_soil[Ssnow_vec>0] = 0 
        ET_vec_soil = ET_vec
        ET_vec_soil[Ssnow_vec>0] = 0
        
        for(i in c(2:length(Ssoil_vec))){
          #calculate soil storage
          if(is.na(Ssoil_vec[i-1])){
            Ssoil_vec[i-1] = 0
          }
          Ssoil_vec[i] = soilstrg_func(Ssoil_vec[i-1], P_vec_soil[i], Ssmax, ET_vec_soil[i])
          soil_overflw[i] = soiloverflow_func(Ssoil_vec[i-1], P_vec[i], Ssmax)
        }
      }
    }
  #all P below critical temperature set to 0, as it falls as snow
  P_vec[Ta_vec<Tcrit_used] = 0
  out_df = data.frame(date_vec, P_vec, Ta_vec, ET_vec, soil_overflw, Ssoil_vec/Ssmax, Ssnow_vec, P_snow_vec)  
  return(out_df)
  },P_vec = as.data.frame(rain_mat), Ta_vec = as.data.frame(Temp_mat), ET_vec = as.data.frame(ET_mat), Ssmax = AWC_in, Tcrit = Tcrit_in, SIMPLIFY = F)
  #change file name according to parameter set?
  #save(class_input_df_list, file = paste0(path_out, "class_input_df_list.Rdata"))
  return(class_input_df_list)
  
} #end of function  
  


#
#  correlation mechanism function -----
#

cor_mech_func = function(){

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Create one dataframe for each catchment with all information -----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  #transfer day of year to a julian date
  doY_to_date = function(in_mat){
    year_dates = as.Date(paste0(in_mat$Year, "-01-01"), format = "%Y-%m-%d")
    in_mat[,"date_Pdaily"] = as.Date(as.numeric(in_mat[,"doY_Pdaily"]), origin = year_dates)-1
    in_mat[,"date_Pweekly"] = as.Date(as.numeric(in_mat[,"doY_Pweekly"]), origin = year_dates)-1
    in_mat[,"date_Peff"] = as.Date(as.numeric(in_mat[,"doY_Peff"]), origin = year_dates)-1
    in_mat[,"date_Psnow"] = as.Date(as.numeric(in_mat[,"doY_Psnow"]), origin = year_dates)-1
    in_mat[,"date_AMAX"] = as.Date(as.numeric(in_mat[,"doY_AMAX"]), origin = year_dates)-1
    return(in_mat)
  }
  
  
  #merge the dataframe by Years so that only years with AMAX available are used for climate calculations
  par_flood_df = mapply(function(amax, rain, week, eff, snow){
    cat_name_AMAX = amax[[1]]
    AMAX_df = amax[[2]]
    years_AMAX = years(as.Date(AMAX_df[,1]))
    date_origin_AMAX = as.Date(paste0(years_AMAX, "-01-01"), format = "%Y-%m-%d")
    AMAX_date = as.Date(as.numeric(AMAX_df[,3]), origin = date_origin_AMAX)-1
    AMAX_df = data.frame(years_AMAX, AMAX_date, as.numeric(AMAX_df[,3]), AMAX_df[,2])
    colnames(AMAX_df) = c("Year", "AMAX_date","doY_AMAX", "AMAX")
    
    rain_df = as.data.frame(rain)
    years_climate = rain_df[,1]
    date_origin = as.Date(paste0(years_climate, "-01-01"), format = "%Y-%m-%d")
    rain_date = as.Date(as.numeric(rain_df$doY), origin = date_origin)-1
    rain_df = data.frame(rain_date, rain_df$doY, rain_df$MaxPdaily)
    
    week_df = as.data.frame(week)
    week_date = as.Date(as.numeric(week_df$doY), origin = date_origin)-1
    week_df = data.frame(week_date, week_df$doY, week_df$MaxPweek)
    
    eff_df = as.data.frame(eff)
    eff_date = as.Date(as.numeric(eff_df$doY), origin = date_origin)-1
    eff_df = data.frame(eff_date, eff_df$doY, eff_df$MaxPeff)
    
    snow_df = as.data.frame(snow)
    snow_date = as.Date(as.numeric(snow_df$doY), origin = date_origin)-1
    snow_df = data.frame(snow_date, snow_df$doY, snow_df$MaxPsnow)
    
    climate_df = data.frame(years_climate, rain_df, week_df, eff_df, snow_df)
    colnames(climate_df) = c("Year", "date_Pdaily",  "doY_Pdaily", "Pdaily", "date_Pweekly", "doY_Pweekly", "Pweekly", "date_Peff", "doY_Peff", "Peff", "date_Psnow", "doY_Psnow", "Psnow")
    
    out = merge(climate_df, AMAX_df, by = "Year", all.x = T)
    list(cat_name_AMAX, out)
    
  }, amax = AMAX_list[cat_ind_rev], rain = P_dailymax, week = P_weeklymax, eff = P_effmax, snow = P_snowmax, SIMPLIFY = F)
  
  save(par_flood_df, file = paste0(file_output_path, "par_flood_df.Rdata"))
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Calculations - Cor month -----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  return_max = function(x){
    isempty_c = x
    isempty_c[is.na(isempty_c)] = 0
    if(all(isempty_c==0)){
      max_val = NA
      #max_date = NA
    } else {
      max_val = max(x)
    }
    max_val
  }
  
  cor_mat = matrix(NA, nrow = length(par_flood_df), ncol = 4)
  colnames(cor_mat) = c("Pdaily", "Pweekly", "Peff", "Psnow")
  
  
  for(i in c(1:length(par_flood_df))){
    count_i = i
    #get AMAX flood dates out of par_flood_df
    temp = par_flood_df[[count_i]][[2]]
    #turn factors to numeric
    ind_col <- sapply(temp, is.factor)
    temp[ind_col] <- lapply(temp[ind_col], function(y) as.numeric(as.character(y)))
    
    AMAX_vals = temp["AMAX"]
    date_end = temp["AMAX_date"]
    
    #temp_out: for every peak flow date, find the maximum mechanism value in the previous time frame 
    temp_out = apply(date_end, 1, function(y){
      #if the date is NA, set output line to NA
      if(is.na(y)){
        temp_vec = rep(NA, 4)
      } else {
        #create sequence of dates based on peak flow date and the threshold for each mechanism
        date_seq = seq(as.Date(y, format = "%Y-%m-%d")-day_thresh, as.Date(y, format = "%Y-%m-%d"), by = "days")
        date_ind = ind_df[,1] %in% date_seq
        
        date_seq_P = seq(as.Date(y, format = "%Y-%m-%d")-day_thresh_Pdaily, as.Date(y, format = "%Y-%m-%d"), by = "days")
        date_ind_P = ind_df[,1] %in% date_seq_P
        
        #create dataframe with daily values selected according to date sequence
        mech_df = data.frame(P_weeklymat[date_ind,count_i], P_effmat[date_ind,count_i], P_snowmat[date_ind,count_i])
        
        #find maximum value in the dataframe (seperate for Pdaily)
        temp_vec1 = return_max(P_dailymat[date_ind_P,count_i])
        temp_vec2 = apply(mech_df,2,return_max)
        temp_vec = c(temp_vec1, temp_vec2)
        names(temp_vec) = c("Pdaily", "Pweekly", "Peff", "Psnow")
        
      }
      return(temp_vec)
    })
    t_temp_out = t(temp_out)
    
    out = apply(t_temp_out, 2, function(x){
      #count how many peak values have form pairs with each mechanism - if not enough values for correlation available, set spearman rank (cor_out) to NA
      pair_count_df = data.frame(x, AMAX_vals)
      pair_count = sum(apply(pair_count_df,1,function(x){all(!is.na(x))}))
      if(pair_count>=cor_thresh){ #only calculate if enough pairs for correlation are available
        cor_out = cor(x, AMAX_vals, "na.or.complete", method = "spearman")
        cor_out = round(cor_out, 2)
      } else {
        cor_out = NA
      }
      cor_out
    })
    cor_mat[count_i,] = out
  }
  
  #set values smaller than threshold (0) to NA
  cor_mat[cor_mat<spearman_thresh] = NA
  
  #find max spearman value/values
  cor_max_list = apply(cor_mat, 1, function(x){
    max_val = max(x, na.rm = T)
    #equal spearman is defined as beeing max or within 0.05 of max value
    max_ind = which(x >= max_val-0.05) #==max_val
    x[max_ind]
  })
  
  cor_mech = lapply(cor_max_list, FUN = function(x){
    length_val = length(x)
    if(length_val==1){
      out = names(x)
      spearman_val = unname(x)
    } else if(length_val>1){
      out = "mixed"
      spearman_val = max(unname(x)) #before: "= NA"
    } else {
      out = "no_class"
      spearman_val = NA
    }
    c(out, spearman_val)
  })
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Merge results to data frame -----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  cat_names = lapply(par_flood_df, function(x){x[[1]]})
  cor_merge_df = data.frame(unlist(cat_names), do.call(rbind, cor_mech), stringsAsFactors = F)
  colnames(cor_merge_df) = c("gsim.no", "cor_mech_month", "spearman_rank")
  cor_merge_df[,3] = as.numeric(cor_merge_df[,3])
  
  write.csv(cor_merge_df, file = paste0(file_output_path, "mech_out.csv"))
  save(cor_merge_df, file = paste0(file_output_path, "mech_out.Rdata"))
  cor_stats = table(cor_merge_df$cor_mech_month)
  save(cor_stats, file = paste0(file_output_path, "stats.Rdata"))
}



#
#  event based classification algorithm  -----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#output is data frame with columns: 


event_classification = function(cat_name,day_thresh, sat_thresh, ts_type = "mean"){
  cat = which(pfdf_names == cat_name) #test for missing flood dates (e.g. are they NA or not there?)
  flood_date = par_flood_df[[cat]][[2]][,"AMAX_date"]
  P_val = P_dailymat[,cat]
  snow_val = P_snowmat[,cat]
  sat_val = P_satmat[,cat]
  
  df = data.frame(ind_df[,1], P_val, snow_val, sat_val)
  colnames(df) = c("date", "P", "S", "sat")

  if(ts_type == "mean"){
    P_thresh = mean(P_val[which(P_val>0)], na.rm = T)
    snow_thresh = mean(snow_val[which(snow_val>0)], na.rm = T)
  } else if(ts_type == "median"){
    P_thresh = median(P_val[which(P_val>0)], na.rm = T)
    snow_thresh = median(snow_val[which(snow_val>0)], na.rm = T)
  } else {
    stop("Pick a threshold type")
  }

  
  # P_thresh = mean(P_val[which(P_val>0)], na.rm = T)
  # snow_thresh = mean(snow_val[which(snow_val>0)], na.rm = T)
  
  #create data frame with relevant information to reach decision of flood mechanism
  decision_df_list = lapply(flood_date, FUN = function(x){
    if(is.na(x)){ #flood date missing
      out = rep(list(NA), 6)
    } else if(!(x %in% ind_df[,1])){ #dates missing
      out = rep(list(NA), 6)
    } else {
      flood_day = which(ind_df[,1] == x)
      temp_df = df[max((flood_day-(day_thresh-1)),1):flood_day,] #max function so that index does not go below 1
      Pmax = max(temp_df$P)
      Smax = max(temp_df$S)
      Pmax_day = temp_df$date[which.max(temp_df$P)]
      #what is soil saturation on day BEFORE max. rainfall
      satPmax = df$sat[max((which(df$date == (Pmax_day-1))),1)] #if flood date is first day of time series use same day saturation
      Smax_day = temp_df$date[which.max(temp_df$S)]
      #How much rainfall on day of max snowmelt
      PSmax = temp_df$P[which(temp_df$date == Smax_day)]
      out = list(round(Pmax,2), Pmax_day, round(satPmax,2), round(Smax,2), Smax_day, round(PSmax,2))
      #if no maximum is found, output is length zero, change to NA
      out = lapply(out, function(x){ifelse(length(x)==0, NA, x)})
    }
    out
  })
  
  decision_df = as.data.frame(rbindlist(decision_df_list))
  colnames(decision_df) = c("Pmax", "Pmax_day", "satPmax", "Smax", "Smax_day", "PSmax")
  
  ELSE = TRUE
  decision_result = decision_df %>% mutate(.,mech_out = with(.,case_when(
    (Smax>snow_thresh & PSmax>P_thresh) ~ "rainsnow",
    (Smax>snow_thresh & PSmax<P_thresh) ~ "snowmelt",
    ((is.na(Smax) | Smax<snow_thresh) & satPmax>sat_thresh & Pmax>P_thresh) ~ "soilsat",
    ((is.na(Smax) | Smax<snow_thresh) & satPmax>sat_thresh & Pmax<P_thresh) ~ "noclass",
    ((is.na(Smax) | Smax<snow_thresh) & satPmax<sat_thresh & Pmax>P_thresh) ~ "rainfall",
    ((is.na(Smax) | Smax<snow_thresh) & satPmax<sat_thresh & Pmax<P_thresh) ~ "noclass",
    ELSE ~ "missingData"
  )))
  
  #turn numbers back into dates
  decision_result = decision_result %>% mutate(Pmax_day = as.Date(Pmax_day,  origin = "1970-01-01"))
  decision_result = decision_result %>% mutate(Smax_day = as.Date(Smax_day,  origin = "1970-01-01"))
  return(decision_result)
}

event_classification_thresholds = function(cat_name,day_thresh, sat_thresh, ts_type = "mean", threshold_min){
  cat = which(pfdf_names == cat_name) #test for missing flood dates (e.g. are they NA or not there?)
  flood_date = par_flood_df[[cat]][[2]][,"AMAX_date"]
  P_val = P_dailymat[,cat]
  snow_val = P_snowmat[,cat]
  sat_val = P_satmat[,cat]
  
  df = data.frame(ind_df[,1], P_val, snow_val, sat_val)
  colnames(df) = c("date", "P", "S", "sat")
  
  if(ts_type == "mean"){
    P_thresh = max(threshold_min, mean(P_val[which(P_val>0)], na.rm = T))
    snow_thresh = max(threshold_min, mean(snow_val[which(snow_val>0)], na.rm = T))
  } else if(ts_type == "median"){
    P_thresh = max(threshold_min, median(P_val[which(P_val>0)], na.rm = T))
    snow_thresh = max(threshold_min, median(snow_val[which(snow_val>0)], na.rm = T))
  } else {
    stop("Pick a threshold type")
  }
  
  
  # P_thresh = mean(P_val[which(P_val>0)], na.rm = T)
  # snow_thresh = mean(snow_val[which(snow_val>0)], na.rm = T)
  
  #create data frame with relevant information to reach decision of flood mechanism
  decision_df_list = lapply(flood_date, FUN = function(x){
    if(is.na(x)){ #flood date missing
      out = rep(list(NA), 6)
    } else if(!(x %in% ind_df[,1])){ #dates missing
      out = rep(list(NA), 6)
    } else {
      flood_day = which(ind_df[,1] == x)
      temp_df = df[max((flood_day-(day_thresh-1)),1):flood_day,] #max function so that index does not go below 1
      Pmax = max(temp_df$P)
      Smax = max(temp_df$S)
      Pmax_day = temp_df$date[which.max(temp_df$P)]
      #what is soil saturation on day BEFORE max. rainfall
      satPmax = df$sat[max((which(df$date == (Pmax_day-1))),1)] #if flood date is first day of time series use same day saturation
      Smax_day = temp_df$date[which.max(temp_df$S)]
      #How much rainfall on day of max snowmelt
      PSmax = temp_df$P[which(temp_df$date == Smax_day)]
      out = list(round(Pmax,2), Pmax_day, round(satPmax,2), round(Smax,2), Smax_day, round(PSmax,2))
      #if no maximum is found, output is length zero, change to NA
      out = lapply(out, function(x){ifelse(length(x)==0, NA, x)})
    }
    out
  })
  
  decision_df = as.data.frame(rbindlist(decision_df_list))
  colnames(decision_df) = c("Pmax", "Pmax_day", "satPmax", "Smax", "Smax_day", "PSmax")
  
  ELSE = TRUE
  decision_result = decision_df %>% mutate(.,mech_out = with(.,case_when(
    (Smax>snow_thresh & PSmax>P_thresh) ~ "rainsnow",
    (Smax>snow_thresh & PSmax<P_thresh) ~ "snowmelt",
    ((is.na(Smax) | Smax<snow_thresh) & satPmax>sat_thresh & Pmax>P_thresh) ~ "soilsat",
    ((is.na(Smax) | Smax<snow_thresh) & satPmax>sat_thresh & Pmax<P_thresh) ~ "noclass",
    ((is.na(Smax) | Smax<snow_thresh) & satPmax<sat_thresh & Pmax>P_thresh) ~ "rainfall",
    ((is.na(Smax) | Smax<snow_thresh) & satPmax<sat_thresh & Pmax<P_thresh) ~ "noclass",
    ELSE ~ "missingData"
  )))
  
  #turn numbers back into dates
  decision_result = decision_result %>% mutate(Pmax_day = as.Date(Pmax_day,  origin = "1970-01-01"))
  decision_result = decision_result %>% mutate(Smax_day = as.Date(Smax_day,  origin = "1970-01-01"))
  return(decision_result)
}

event_classification_parts = function(cat_name,day_thresh, sat_thresh, part_ts = 2/3){
  cat = which(pfdf_names == cat_name) 
  flood_date = par_flood_df[[cat]][[2]][,"AMAX_date"]
  P_val = P_dailymat[,cat]
  snow_val = P_snowmat[,cat]
  sat_val = P_satmat[,cat]
  
  df = data.frame(ind_df[,1], P_val, snow_val, sat_val)
  colnames(df) = c("date", "P", "S", "sat")
  
  P7_mean = mean(rollapply(P_val, day_thresh, sum, partial = T, align = "right"))
  
  #create data frame with relevant information to reach decision of flood mechanism
  decision_df_list = lapply(flood_date, FUN = function(x){
    if(is.na(x)){ #flood date missing
      out = rep(list(NA), 9)
    } else if(!(x %in% ind_df[,1])){ #dates missing
      out = rep(list(NA), 9)
    } else {
      flood_day = which(ind_df[,1] == x)
      temp_df = df[max((flood_day-(day_thresh-1)),1):flood_day,] #max function so that index does not go below 1
      Pmax = max(temp_df$P)
      P7 = sum(temp_df$P)
      S7 = sum(temp_df$S)
      Ptotal = sum(c(P7, S7), na.rm = T)
      Pmax_day = temp_df$date[which.max(temp_df$P)]
      #what is soil saturation on day BEFORE max. rainfall
      satPmax = df$sat[max((which(df$date == (Pmax_day-1))),1)] #if flood date is first day of time series use same day saturation
      satQmax = temp_df$sat[day_thresh] 
      Smax_day = temp_df$date[which.max(temp_df$S)]
      #How much rainfall on day of max snowmelt
      #PSmax = temp_df$P[which(temp_df$date == Smax_day)]
      out = list(round(P7_mean), round(Ptotal, 2), round(P7, 2), round(Pmax,2), Pmax_day, round(satPmax,2), round(satQmax,2), round(S7,2), Smax_day)
      #if no maximum is found, output is length zero, change to NA
      out = lapply(out, function(x){ifelse(length(x)==0, NA, x)})
    }
    out
  })
  
  decision_df = as.data.frame(rbindlist(decision_df_list))
  colnames(decision_df) = c("P7_mean", "Ptotal", "P7", "Pmax", "Pmax_day", "satPmax", "satQmax", "S7", "Smax_day")
  
  ELSE = TRUE
  decision_result = decision_df %>% mutate(.,mech_out = with(.,case_when(
    #minimum rainfall limit
    (Ptotal < 1) ~ "noclass",
    (P7>=(part_ts*Ptotal) & satPmax>=sat_thresh & Pmax>=(part_ts*P7)) ~ "rainfall",
    (P7>=(part_ts*Ptotal) & satPmax>=sat_thresh & Pmax<(part_ts*P7) & P7>=P7_mean) ~ "soilsat",
    (P7>=(part_ts*Ptotal) & satPmax>=sat_thresh & Pmax<(part_ts*P7) & P7<P7_mean) ~ "noclass",
    (P7>=(part_ts*Ptotal) & satPmax<sat_thresh & Pmax>=(part_ts*P7)) ~ "rainfall",
    (P7>=(part_ts*Ptotal) & satPmax<sat_thresh & Pmax<(part_ts*P7) & satQmax>=sat_thresh & P7>=P7_mean) ~ "soilsat",
    (P7>=(part_ts*Ptotal) & satPmax<sat_thresh & Pmax<(part_ts*P7) & satQmax>=sat_thresh & P7<P7_mean) ~ "noclass",
    (P7>=(part_ts*Ptotal) & satPmax<sat_thresh & Pmax<(part_ts*P7) & satQmax<sat_thresh) ~ "noclass",
    (P7<(part_ts*Ptotal) & S7>=(part_ts*Ptotal)) ~ "snowmelt",
    (P7<(part_ts*Ptotal) & S7<(part_ts*Ptotal)) ~ "rainandsnow",
    ELSE ~ "missingData"
  )))
  
  #turn numbers back into dates
  decision_result = decision_result %>% mutate(Pmax_day = as.Date(Pmax_day,  origin = "1970-01-01"))
  decision_result = decision_result %>% mutate(Smax_day = as.Date(Smax_day,  origin = "1970-01-01"))
  return(decision_result)
}

event_classification_test = function(cat_name,day_thresh, sat_thresh, part_ts = 2/3){
  cat = which(pfdf_names == cat_name) 
  flood_date = par_flood_df[[cat]][[2]][,"AMAX_date"]
  magn = par_flood_df[[cat]][[2]][,"AMAX"]
  P_val = P_dailymat[,cat]
  snow_val = P_snowmat[,cat]
  sat_val = P_satmat[,cat]
  
  
  df = data.frame(ind_df[,1], P_val, snow_val, sat_val)
  colnames(df) = c("date", "P", "S", "sat")
  
  P7_mean = mean(rollapply(P_val, day_thresh, sum, partial = T, align = "right"))
  P_wetdays = P_val[P_val>1]
  Pextreme_ts = rev(sort(P_wetdays))[round(length(P_wetdays)*0.01)]
  
  #create data frame with relevant information to reach decision of flood mechanism
  decision_df_list = lapply(flood_date, FUN = function(x){
    if(is.na(x)){ #flood date missing
      out = rep(list(NA), 9)
    } else if(!(x %in% ind_df[,1])){ #dates missing
      out = rep(list(NA), 9)
    } else {
      flood_day = which(ind_df[,1] == x)
      temp_df = df[max((flood_day-(day_thresh-1)),1):flood_day,] #max function so that index does not go below 1
      Pmax = max(temp_df$P)
      P7 = sum(temp_df$P)
      S7 = sum(temp_df$S)
      Ptotal = sum(c(P7, S7), na.rm = T)
      Pmax_day = temp_df$date[which.max(temp_df$P)]
      #what is soil saturation on day BEFORE max. rainfall
      satPmax = df$sat[max((which(df$date == (Pmax_day-1))),1)] #if flood date is first day of time series use same day saturation
      satQmax = temp_df$sat[day_thresh] 
      Smax_day = temp_df$date[which.max(temp_df$S)]
      #How much rainfall on day of max snowmelt
      #PSmax = temp_df$P[which(temp_df$date == Smax_day)]
      out = list(round(Pextreme_ts),round(P7_mean), round(Ptotal, 2), round(P7, 2), round(Pmax,2), Pmax_day, round(satPmax,2), round(satQmax,2), round(S7,2), Smax_day)
      #if no maximum is found, output is length zero, change to NA
      out = lapply(out, function(x){ifelse(length(x)==0, NA, x)})
    }
    out
  })
  
  decision_df = as.data.frame(rbindlist(decision_df_list))
  decision_df = data.frame(flood_date, magn, decision_df)
  colnames(decision_df) = c("flood_date", "AMAX", "Pextreme", "P7_mean", "Ptotal", "P7", "Pmax", "Pmax_day", "satPmax", "satQmax", "S7", "Smax_day")
  #colnames(decision_df) = c("Pextreme", "P7_mean", "Ptotal", "P7", "Pmax", "Pmax_day", "satPmax", "satQmax", "S7", "Smax_day")
  
  ELSE = TRUE
  decision_result = decision_df %>% mutate(.,mech_out = with(.,case_when(
    #minimum rainfall limit
    (Ptotal < 1) ~ "noclass",
    (P7>=(part_ts*Ptotal) & satPmax>=sat_thresh & Pmax>=(part_ts*P7)) ~ "rainfall",
    (P7>=(part_ts*Ptotal) & satPmax>=sat_thresh & Pmax<(part_ts*P7) & P7>=P7_mean) ~ "soilsat",
    (P7>=(part_ts*Ptotal) & satPmax>=sat_thresh & Pmax<(part_ts*P7) & P7<P7_mean) ~ "noclass",
    (P7>=(part_ts*Ptotal) & satPmax<sat_thresh & Pmax>=(part_ts*P7)) ~ "rainfall",
    (P7>=(part_ts*Ptotal) & satPmax<sat_thresh & Pmax<(part_ts*P7) & satQmax>=sat_thresh & P7>=P7_mean) ~ "soilsat",
    (P7>=(part_ts*Ptotal) & satPmax<sat_thresh & Pmax<(part_ts*P7) & satQmax>=sat_thresh & P7<P7_mean) ~ "noclass",
    (P7>=(part_ts*Ptotal) & satPmax<sat_thresh & Pmax<(part_ts*P7) & satQmax<sat_thresh) ~ "noclass",
    (P7<(part_ts*Ptotal) & S7>=(part_ts*Ptotal)) ~ "snowmelt",
    (P7<(part_ts*Ptotal) & S7<(part_ts*Ptotal)) ~ "rainandsnow",
    ELSE ~ "missingData"
  )))
  
  #turn numbers back into dates
  decision_result = decision_result %>% mutate(Pmax_day = as.Date(Pmax_day,  origin = "1970-01-01"))
  decision_result = decision_result %>% mutate(Smax_day = as.Date(Smax_day,  origin = "1970-01-01"))
  return(decision_result)
}

event_classification_freudiger = function(cat_name,day_thresh, sat_thresh){
  cat = which(pfdf_names == cat_name) #test for missing flood dates (e.g. are they NA or not there?)
  flood_date = par_flood_df[[cat]][[2]][,"AMAX_date"]
  P_val = P_dailymat[,cat]
  snow_val = P_snowmat[,cat]
  sat_val = P_satmat[,cat]
  strg_val = P_snowstrgmat[,cat]
  
  df = data.frame(ind_df[,1], P_val, snow_val, sat_val, strg_val)
  colnames(df) = c("date", "P", "S", "sat", "storage")
  # 
  # if(ts_type == "mean"){
  #   P_thresh = max(5, mean(P_val[which(P_val>0)], na.rm = T))
  #   snow_thresh = max(5, mean(snow_val[which(snow_val>0)], na.rm = T))
  # } else if(ts_type == "median"){
  #   P_thresh = max(5, median(P_val[which(P_val>0)], na.rm = T))
  #   snow_thresh = max(5, median(snow_val[which(snow_val>0)], na.rm = T))
  # } else {
  #   stop("Pick a threshold type")
  # }
  # 
  
  
  P_thresh = max(5, mean(P_val[which(P_val>0)], na.rm = T))
  RoS_thresh = 3
  strg_thresh = 10
  snow_percent = 0.2
  snow_thresh = max(5, mean(snow_val[which(snow_val>0)], na.rm = T))
  
  #create data frame with relevant information to reach decision of flood mechanism
  decision_df_list = lapply(flood_date, FUN = function(x){
    if(is.na(x)){ #flood date missing
      out = rep(list(NA), 7)
    } else if(!(x %in% ind_df[,1])){ #dates missing
      out = rep(list(NA), 7)
    } else {
      flood_day = which(ind_df[,1] == x)
      temp_df = df[max((flood_day-(day_thresh-1)),1):flood_day,] #max function so that index does not go below 1
      Pmax = max(temp_df$P)
      Smax = max(temp_df$S)
      Pmax_day = temp_df$date[which.max(temp_df$P)]
      #what is soil saturation on day BEFORE max. rainfall
      satPmax = df$sat[max((which(df$date == (Pmax_day-1))),1)] #if flood date is first day of time series use same day saturation
      Smax_day = temp_df$date[which.max(temp_df$S)]
      #How much rainfall on day of max snowmelt
      PSmax = temp_df$P[which(temp_df$date == Smax_day)]
      strg_snow = df$storage[max((which(df$date == (Smax_day-1))),1)]
      #How much storage on day BEFORE max snowmelt
      out = list(round(Pmax,2), Pmax_day, round(satPmax,2), round(Smax,2), Smax_day, round(PSmax,2), round(strg_snow, 2))
      #if no maximum is found, output is length zero, change to NA
      out = lapply(out, function(x){ifelse(length(x)==0, NA, x)})
    }
    out
  })
  
  decision_df = as.data.frame(rbindlist(decision_df_list))
  colnames(decision_df) = c("Pmax", "Pmax_day", "satPmax", "Smax", "Smax_day", "PSmax", "strg_snow")
  decision_df = decision_df %>% mutate(
    perc_snow = Smax/(PSmax+0.0001) #to avoid division by zero
  )
  decision_df[!is.finite(decision_df$perc_snow), "perc_snow"] = NA #turn NaN and Inf into NA
  
  ELSE = TRUE
  decision_result = decision_df %>% mutate(.,mech_out = with(.,case_when(
    (strg_snow>strg_thresh & PSmax>RoS_thresh & perc_snow>=snow_percent) ~ "rainsnow",
    (Smax>snow_thresh & PSmax<RoS_thresh) ~ "snowmelt",
    ((is.na(Smax) | Smax<snow_thresh) & satPmax>sat_thresh & Pmax>P_thresh) ~ "soilsat",
    ((is.na(Smax) | Smax<snow_thresh) & satPmax>sat_thresh & Pmax<P_thresh) ~ "noclass",
    ((is.na(Smax) | Smax<snow_thresh) & satPmax<sat_thresh & Pmax>P_thresh) ~ "rainfall",
    ((is.na(Smax) | Smax<snow_thresh) & satPmax<sat_thresh & Pmax<P_thresh) ~ "noclass",
    ELSE ~ "missingData"
  )))
  
  #turn numbers back into dates
  decision_result = decision_result %>% mutate(Pmax_day = as.Date(Pmax_day,  origin = "1970-01-01"))
  decision_result = decision_result %>% mutate(Smax_day = as.Date(Smax_day,  origin = "1970-01-01"))
  return(decision_result)
}


event_classification_collins = function(cat_name,day_thresh=2, sat_thresh, ts_type = "mean"){
  cat = which(pfdf_names == cat_name) #test for missing flood dates (e.g. are they NA or not there?)
  flood_date = par_flood_df[[cat]][[2]][,"AMAX_date"]
  P_val = P_dailymat[,cat]
  snow_val = P_snowmat[,cat]
  sat_val = t(soil_sat)[,cat]
  
  df = data.frame(ind_df[,1], P_val, snow_val, sat_val)
  colnames(df) = c("date", "P", "S", "sat")
  
  if(ts_type == "mean"){
    P_thresh = mean(P_val[which(P_val>0)], na.rm = T)
    snow_thresh = mean(snow_val[which(snow_val>0)], na.rm = T)
  } else if(ts_type == "median"){
    P_thresh = median(P_val[which(P_val>0)], na.rm = T)
    snow_thresh = median(snow_val[which(snow_val>0)], na.rm = T)
  } else {
    stop("Pick a threshold type")
  }
  
  
  # P_thresh = mean(P_val[which(P_val>0)], na.rm = T)
  # snow_thresh = mean(snow_val[which(snow_val>0)], na.rm = T)
  
  #create data frame with relevant information to reach decision of flood mechanism
  decision_df_list = lapply(flood_date, FUN = function(x){
    if(is.na(x)){ #flood date missing
      out = rep(list(NA), 3)
    } else if(!(x %in% ind_df[,1])){ #dates missing
      out = rep(list(NA), 3)
    } else {
      flood_day = which(ind_df[,1] == x)
      temp_df = df[(flood_day-1):(flood_day+1),] #max function so that index does not go below 1
      P3d = sum(temp_df$P)
      S3d = sum(temp_df$S)
      temp_df2 = df[max((flood_day-(13)),1):(flood_day+1),]
      S14d = sum(temp_df2$S)
      out = list(round(P3d,2), round(S3d,2), round(S14d,2))
      #if no maximum is found, output is length zero, change to NA
      #out = lapply(out, function(x){ifelse(length(x)==0, NA, x)})
    }
    out
  })
  
  decision_df = as.data.frame(rbindlist(decision_df_list))
  colnames(decision_df) = c("Pval", "Sval", "Sval_long")
  
  ELSE = TRUE
  decision_result = decision_df %>% mutate(.,mech_out = with(.,case_when(
    (Sval>0 & Pval <1) ~ "snowmelt",
    (Sval>0 & Pval>0 & Pval<(2*Sval)) ~ "rainsnow",
    (Sval_long>0 & Pval>0 & Pval<(2*Sval_long)) ~ "rain/snowmelt",
    (Pval>(2*Sval_long)) ~ "rain",
    ELSE ~ "missingData"
  )))
  return(decision_result)
}


# Event classification optimised for sensitivity analysis ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

event_classification_df = function(cat_name, day_thresh, quant_thresh_rain, quant_thresh_snow, df_list, full_cat_list = pfdf_names, flood_df_in = flood_df, date_vec){
  #date_vec: time stamp for timeseries
  cat = which(full_cat_list == cat_name) 
  flood_date = .subset2(flood_df_in[[cat]][[2]], "AMAX_date")
  flood_magn = as.numeric(as.character(.subset2(flood_df_in[[cat]][[2]], "AMAX")))
  
  flood_data_df = data.frame(flood_date, flood_magn)
  colnames(flood_data_df) = c("flood_date", "flood_magn")
  
  df = df_list[[cat]]
  P_val = .subset2(df, "P_vec")
  snow_val = .subset2(df, "P_snow_vec")
  sat_val = .subset2(df, "Ssoil_vec.Ssmax")
  dates_full = .subset2(df, 1)
  
  #due to the structure of the classification function all decision thresholds need to be the same length as the data, therefore repeated 34 times
  P7_rollsum = RcppRoll::roll_sum(c(numeric(day_thresh-1), P_val), n = day_thresh, align = "right")
  #add numeric zeros to achieve partial rollsums at start of ts
  P7_mean = mean(P7_rollsum, na.rm = T)
  P7_quant = quantile(P7_rollsum, quant_thresh_rain, na.rm = T)
  Pday_quant = quantile(P_val, quant_thresh_rain, na.rm = T)
  S7_rollsum = RcppRoll::roll_sum(c(numeric(day_thresh-1), snow_val), n = day_thresh, align = "right")
  S7_quant = quantile(S7_rollsum[S7_rollsum>1], quant_thresh_snow, na.rm = T)
  
  temp = round(c(Pday_quant, P7_mean, P7_quant, S7_quant), 2)
  temp_out = matrix(rep(temp, length(flood_magn)), nrow = length(flood_magn), byrow = T) #changed '34' to 'length(flood_magn)' 20.08.19 
  colnames(temp_out) = c("Pday_quant", "P7_mean", "P7_quant", "S7_quant")
  
  
  #create data frame with relevant information to reach decision of flood mechanism
  decision_df_list = lapply(flood_date, FUN = function(x){
    if(is.na(x)){ #flood date missing
      out = rep(list(NA), 6)
    } else if(all(is.na(sat_val))){ #max soil routine is NA
      out = rep(list(NA), 6)
    } else if(!(x %in% date_vec)){ #dates missing
      out = rep(list(NA), 6)
    } else{
      flood_day = which(dates_full == x)
      flood_period = c(1:length(dates_full)) %in% c(max((flood_day-(day_thresh-1)),1):flood_day) #max function so that index does not go below 1
      
      #temp_df = df[max((flood_day-(day_thresh-1)),1):flood_day,] 
      Pmax = max(P_val[flood_period])
      P7 = sum(P_val[flood_period])
      S7 = sum(snow_val[flood_period])
      Ptotal = sum(c(P7, S7), na.rm = T)
      Pmax_frac = Pmax/P7
      sat_start = sat_val[flood_period][1] #saturation on first day of period
      out = list(sat_start, Pmax, P7, Pmax_frac, S7, Ptotal)
      out = lapply(out, round, 2)
      #if no maximum is found, output is length zero, change to NA
      out = lapply(out, function(x){ifelse(length(x)==0, NA, x)})
    }
    out
  })
  
  decision_df_temp = as.data.frame(rbindlist(decision_df_list))
  colnames(decision_df_temp) = c("sat_start", "Pmax", "P7", "Pmax_frac", "S7", "Ptotal")
  decision_df = data.frame(flood_data_df, decision_df_temp, temp_out)
  decision_df
}


event_classification_quantile = function(decision_df, parts_rainsnow = 1/3, parts_fracextreme = 2/3, sat_thresh = 0.9){
  ELSE = TRUE
  decision_result = decision_df %>% mutate(.,mech_out = with(.,case_when(
    #minimum rainfall limit
    (Ptotal < 1) ~ "noclass",
    (!is.na(S7) & (S7/Ptotal)>=parts_rainsnow & (P7/Ptotal)>=parts_rainsnow) ~ "rainandsnow",
    (!is.na(S7) & S7>=S7_quant) ~ "snowmelt",
    (sat_start >= sat_thresh & P7 >= P7_mean) ~ "soilsat",
    (P7 >= P7_quant & Pmax_frac >= parts_fracextreme) ~ "rainfall",
    (P7 >= P7_quant & Pmax_frac < parts_fracextreme) ~ "longrainfall",
    (Pmax >= Pday_quant) ~ "rainfall",
    (Pmax < Pday_quant) ~ "noclass",
    ELSE ~ "missingData"
  )))
  
  #turn numbers back into dates
  decision_result = decision_result %>% mutate(flood_date = as.Date(flood_date,  origin = "1970-01-01"))
  return(decision_result)
}




# Dominant mechanism calculation  ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#old dom mech function----

# dom_mech_func = function(event_mech_df){
#   temp = table(event_mech_df[,"mech_out"])
#   #as which mechanisms are most events classified (with a minimum of ten)
#   mech = names(temp)
#   vals = unname(temp)
#   
#   mech_ind = mech %in% c("rainfall", "soilsat", "rainandsnow", "snowmelt", "longrainfall")
#   
#   if(length(mech_ind)==1){ #if all one mechanism
#     dom_mech = mech
#   } else if(length(vals[mech_ind])==0) { #if only mechanisms are missing and noclass
#     max_ind_na = which.max(vals[!mech_ind])
#     dom_mech = mech[!mech_ind][max_ind_na]
#   } else {
#     if(max(vals[mech_ind])<10){ #(one value which is not noclass and not missingData is classified at least ten times)
#       if(sum(vals[mech_ind])<20){
#         #more missing data or not classified years
#         max_ind_na = which.max(vals[!mech_ind])
#         dom_mech = mech[!mech_ind][max_ind_na]
#       } else {
#         #sum of all mechanisms together is at least 20
#         dom_mech = "mixed"
#       }
#     } else {
#       #which mechanism classifies at least ten years and is maximum?
#       max_ind = which.max(vals[mech_ind])
#       dom_mech = mech[mech_ind][max_ind]
#     }
#   }
#   dom_mech
# }
#####

dom_mech_func = function(event_mech_df){
  mech_vals = c("noclass","soilsat",  "rainfall", "longrainfall", "snowmelt", "rainsnow")
  temp =event_mech_df[,"mech_out"]
  temp = factor(temp, levels =c("missingData", "noclass","soilsat",  "rainfall", "longrainfall", "snowmelt", "rainandsnow"))
  freq_table = table(temp)
  if(sum(freq_table[2:7])>=20){
    out = mech_vals[which.max(freq_table[2:7])]
    sorted_vals = sort(freq_table[2:7], decreasing = T)
    certainty = unname(sorted_vals[2]/sorted_vals[1])
  } else {
    out = "missingData"#"notenoughData"
    certainty = NA
  }
  
  data.frame(out, certainty, stringsAsFactors = F)
}



#
#  plotting function -----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

plot_flood_event_ts = function(cat_name){
  cat = which(pfdf_names == cat_name)
  flood_mag = par_flood_df[[cat]][[2]][,"AMAX"]
  flood_mag = as.numeric(as.character(flood_mag))
  flood_date = par_flood_df[[cat]][[2]][,"AMAX_date"]
  flood_ind = ind_df$date %in% flood_date
  
  P_val = P_dailymat[,cat]
  #eff_val = P_effmat[,cat]
  snow_val = P_snowmat[,cat]
  sat_val = t(soil_sat)[,cat]
  temp_val = rep(-8, length(P_val))
  #flood_mag_vec = rep(0, length(P_val))
  #flood_mag_vec[flood_ind] = flood_mag
  
  df1 = data.frame(ind_df$date, P_val, snow_val, sat_val, temp_val)
  colnames(df1) = c("date", "P", "snow", "sat", "temp")
  
  mag_df = data.frame(c(1:length(flood_mag)), flood_mag)

  plot_list_all = lapply(mag_df[,1], function(i){
    p1 = ggplot(df1, aes(x = date, y = P))+
      geom_line(colour = "red", size= 1)+
      geom_line(data = df1, aes(x = date, y = snow), colour = "blue", size= 1)+
      #geom_line(data = df1, aes(x = date, y = long), colour = "orange", size= 1)+
      #geom_line(data = df1, aes(x = date, y = eff), colour = "green", size= 1)+
      geom_line(data = df1, aes(x = date, y = temp, colour = sat), size= 2)+scale_colour_gradient(low = "red", high = "blue")+
      geom_dotplot(data = df1[flood_ind,], aes(x = date), colour = NA , dotsize = 0.5, binwidth = 1)+
      #geom_segment(data = df1[flood_ind,], aes(xend = date, yend = mag))+
      scale_x_date(limits = c(flood_date[i]-7, flood_date[i]+1))+
      labs(title = paste(flood_date[i], flood_mag[i], "m3/s"))+
      ylim(-10, max(c(max(P_val, na.rm = T), max(snow_val, na.rm = T)), na.rm = T))+
      guides(colour=FALSE)
    p1
  })
  
  do.call(grid.arrange, plot_list_all)
}


#
#  multiplot function -----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#source
#http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}



#function to group catchments by 2x2 degree cell ----
#source: https://gis.stackexchange.com/questions/48416/aggregating-points-to-grid-using-r
ji <- function(xy, origin=c(0,0), cellsize=c(2,2)) {
  t(apply(xy, 1, function(z) cellsize/2+origin+cellsize*(floor((z - origin)/cellsize))))
}



