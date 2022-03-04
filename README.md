# DGVMresilience.R
R script for analysing global vegetation resilience trends in a Dynamic General Vegetation Model (DGVM), using STL decomposition for seasonal detrending of input data (default net primary productivity), autocorrelation (or standard deviation) as a proxy for resilience, and the Mann Kendall trend test for trend strength (see Lenton et al., 2022 for details).

Code author: Dr. David I. Armstrong McKay [Georesilience Analytics; Uni. Exeter]

Written for: Lenton et al., 2022: "A resilience sensing system for the biosphere", PhilTransB, https://doi.org/10.1098/rstb.2021.0383 - used to generate Figure 1

ISIMIP2b data available from: https://doi.org/10.5880/PIK.2019.012 and https://data.isimip.org/

Some bugs remain for the parallel processing component - fix will be made available in future updates (and serial processing is stable in the meantime)
