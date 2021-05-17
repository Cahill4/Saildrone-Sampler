# Saildrone-Sampler
The following Python code uses CFSv2 weather data to replicate Saildrone data retrieval. Much of this code is built upon UW's, Shuyi Chen, pyhycom code.

# Data
- All data is from the NOAA/NCEP modified Coupled Forecast System v2 Model
- Currently only ocean data is utilized but that may change
- Data spans two years (2002-2003)
- Data is NetCDF4

# What Code Analyzes
- Goal is to output variables Saildrones may be able to capture and test efficiency of data retireval at different paths
- Output would be valuable as Synthetic Observations for Observing System Simualtion Experiments (OSSEs) 
- Study region is the East Pacific Hurricane Genesis Region (Between 10°N-15°N and 110°W-100°W)
- Two differnt Saildrones are tested at once
- Three different paths can be tested: (Meridional, Zonal, and Cross)
- Options are only available for a 12 day interval but changes can be made for other intervals
- Current code analyzes Saildrone at a constant speed, but updates to come

# Output
- Saildrone trajectories are calcuated from the wind speed and currents at each hourly time-step
- Plot 1: Ocean Temp and Salinty
  - Top two graphs are day one averages overlayed by the paths of each Saildrone
  - Bottomm four are temp and sailinty cross sections over Saildrone's timescale
- Plot 2: Surface and Latent Heat Flux
  - Top two graphs are day one averages overlayed by the paths of each Saildrone
  - Bottom two aretime series plots based on the Saildrone's path
- Plot 3: Wind Stress and Surface Current Speeds (same format as plot 2)
- Plot 4: Precip and Ocean Surface Height (same format as plot 2)

# Acknowledgements
- Much of this code is built off of Shuyi Chen's pyhycom's code - Thank you
- Special thanks to my three incredible mentors - Meghan Cronin, Dongxiao Zhang, and Samantha Wills
- Thank you to entire OCS team - Patrick Berk, Nathan Anderson, Jessica Masich, Annika Margevich, and Ilyana Collins
- Thank you to NCEP for supplying the model data
- Thank you to the NOAA Hollings Program
- NCdump function courtesy of Chris Slocum of Colorado State

