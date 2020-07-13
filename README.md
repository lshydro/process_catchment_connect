# process_catchment_connect
Code and data used analyse the connection between flood processes and catchment characteristics

Calculation_visualisation_process_catchment.r includes all calculations and visualisations to interpret the results. 

function_file_global_flood_classification.R contains the code for the global flood classification as described by Stein et al. (2020)

Data requirements to execute the code: 
- CAMELS dataset by Addor et al. (2017) available at https://ral.ucar.edu/solutions/products/camels
- aws0200NATSGO.csv Available Water Storage data calculated for the CAMELS catchments with data from the NATSGO soil data base. Available as supplemental data for the article submitted to Water Resources Research.


Addor, N., Newman, A.J., Mizukami, N. and Clark, M.P., 2017. The CAMELS data set: catchment attributes and meteorology for large-sample studies. Hydrology and Earth System Sciences (HESS), 21(10), pp.5293-5313.
Stein, L., Pianosi, F. and Woods, R., 2020. Event‐based classification for global study of river flood generating processes. Hydrological Processes, 34(7), pp.1514-1529.
