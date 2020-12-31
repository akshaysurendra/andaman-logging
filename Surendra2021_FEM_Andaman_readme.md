# README FILE

This Surendra2021_FEM_Andaman_readme.txt file was generated on 2020-12-31 by Akshay Surendra

## GENERAL INFORMATION

  1. Title: Data associated with "Varying impacts of logging frequency on tree communities and carbon storage across evergreen and deciduous tropical forests in the Andaman Islands, India" 
  2. Authors: Akshay Surendra (akshaysurendra1@gmail.com), Anand M Osuri, Jayashree Ratnam
  3. Date of data collection: November 2017 to May 2018
  4. Geographic location of data collection: Latitude: 12.13N - 12.40N, Longitude: 92.73E - 92.90E, Middle and North Andaman, Andaman and Nicobar Islands, India
  5. Funding sources that supported data collection: (1) Department of Atomic Energy (DAE), Govt. of India via National Centre for Biological Sciences- Tata Institute of Fundamental Research (NCBS-TIFR) (2) Sir Dorabji Tata Trust

## SHARING/ACCESS INFORMATION

The data can be used for non-commercial purposes (The MIT Licence)
Publication: https://doi.org/10.1016/j.foreco.2020.118791
Data is also hosted on Github along with the code: https://github.com/akshaysurendra/andaman-logging

Recommended citation for this dataset: Surendra A, Osuri M A, Ratnam J (2021) Data from: Varying impacts of logging frequency on tree communities and carbon storage across evergreen and deciduous tropical forests in the Andaman Islands, India. Dryad Digital Repository. https://doi.org/10.5061/dryad.1c59zw3rp

## DATA & FILE OVERVIEW

#### File List:

  1. *plotdata.csv*: This file contains vegetation data corresponding to each plot across all three subplots (subplot A - 50x10m, subplotB - 10x10m and subplotC - 70x70m); subplots A and B were used in the associated manuscript but data from subplotC have not been analyzed
  2. *plot_characteristics.csv*: This file contains slope, elevation, distance to the forest edge, distance to nearest neighbour, temperature (Worldclim) and rainfall (Worldclim) values, and the centroid (lat-long) for each sampled plot
  3. *canopy.csv*: This file contains the canopy cover measured as presence (1) or absence (0) in 45 nearly 1-m intervals along the long-axis of the 50x10m plot; each column corresponds to the total counts in a 5m interval (maximum: 5 and minimum: 0), at two heights: <15m (L) and >20m (H)
  4. *treecodes.csv*: This file contains information on each (morpho-)species encountered in our plots

**PlotID** is a key that can be used to join files A, B and C; file D can be joined with A using the code (and code_new) column

**Additional related data (not included here)**:Total soil carbon % and total soil nitrogen % were measured in each plot at two depths; samples were collected at three equally spaced intervals along the long-axis of the 50x10m plot and pooled before analysis (Ignition-based loss-of-mass estimation); data are available on request.
This is the only version of the dataset used in the analysis associated with the manuscript listed above.

**Size and shape of plots**: Each plot consists of three subplots that are nested within each other and the plots are aligned along a logged stump. In relatively undisturbed baseline plots (B), the plots are oriented in a random direction starting at a large tree of loggable girth. The outermost subplot (plotC) is 70x70m across and encircles a logging gap exactly at its centre; in this plot all trees >=180cm girth are measured and identified as are any stumps that have been logged.

Within this, and 35m from two opposite edges of plotC is plotA, 50x10m subplot that starts exactly at the logging stump and extends for 50m along the central line with 5m on either side; this plot is designed to best capture logging disturbance and all stems >=30cm girth are measured. in the first one-fifth of this 50x10m plot, we placed plotB - a 10x10m subplot within which all stems >=10cm girth were measured. Only data from plotB and plotA were used for the article associated with this paper.

## METHODOLOGICAL INFORMATION

See methods section of https://doi.org/10.1016/j.foreco.2020.118791 ; plotC data has not been included in the manuscript. Plot C is 70x70m (0.49ha) in size and encircles the other two plots right at the centre, and all trees (and cut stumps) above 180cm girth within this larger plot have been have been identified and measured. 

## (1) DATA-SPECIFIC INFORMATION FOR: plotdata.csv
  
  1. Number of variables: 20
  2. Number of cases/rows: 4303
  3. Missing data code: NA
  4. Variable list and description:
      1. treatment - three levels that indicate the logging history experienced by the plot B (baseline) L1 (once-logged) and L2 (twice-logged)
      2. forest_type - deciduous or evergreen forest type
      3. plot_ID - unique ID for each of the plots
      4. gbh_cm - tree girth ~1.3m from the ground (in cms)
      5. gbh_cm_stem2 - tree girth ~1.3m from the ground of the second stem if present (in cms)
      6. gbh_cm_stem3 - tree girth ~1.3m from the ground of the third stem if present (trees with >3 stems are collapsed into this 3rd stem) (in cms)
      7. height_m - height measured in the field using a laser rangefinder - visually estimated for stems <6m (in metres)
      8. code: unique code for each species with first 3 letters indicate Genus and the next three indicate Species in most cases
      9. species_ID - scientific name when known local name if not or a unique morphospecies code
      10. family - Family to which that stem belongs to
      11. ID_level - S: identified to species level G: identified to Genus level F: identified to Family level
      12. plotA - any stem >30cm GBH measured within the 50x10m plot area (adult-trees)
      13. plotB - any stem in the 10-30cm GBH range measured within the 10x10m plot area (juvenile-trees)
      14. plotC - any stem >=180cm GBH measured within the 70x70cm plot area (large-trees unpublished data)
      15. wd_m - wood density using sp_m when available and ge_m if not and if both are absent plot-averaged values (g/cm^3)
      16. wd_sd - variation in the wood density carried forward from species genus or plot-specific values
      17. WD_level - ID level at which wood density was calculated (S: species-specific G: genus-specific P: plot-specific)
      18. agb_Mg - agb / (10^6) - tonnes of carbon
      19. type_plot - B: stem measured in the innermost 10x10m subplot; A: adult tree measured in remaining area of the treeline subplot (40x10m with stems in the first 10x10m already measured in plot B); C: large tree (>=180cm GBH) measured in the outermost 70x70m subplot circumambulating the other two subplots excluding large trees that were already measured in the inner plots
      20. gf - T: tree; P: palm or cane (monocot)

## (2) DATA-SPECIFIC INFORMATION FOR: plot_characteristics.csv
  
  1. Number of variables: 10
  2. Number of cases/rows: 76
  3. Missing data code: NA
  4. Variable list and description:
      1. plotID - unique ID for each of the plots
      2. forest_type - evergreen and deciduous forest type
      3. latitude longitude - as is
      4. elevSTRM - elevation extracted from SRTM satellite data (30m resolution) in meters
      5. slopeSRTM - slope extracted from SRTM satellite data (30m resolution) in degrees
      6. distanceNearestNeighbour_m - distance of the nearest neighbouring plot in meters
      7. LoggingTreatment - three levels that indicate the logging history experienced by the plot B (baseline) L1 (once-logged) and L2 (twice-logged) same as in plotdata.csv
      8. distSettlement_m - distance of the plot to the forest edge in meters
      9. temperature - average annual temperature extracted from WorldClim data (in degree Celsius)
      10. rainfall - total annual rainfall extracted from WorldClim data (in mm)

## (3) DATA-SPECIFIC INFORMATION FOR: canopy.csv

  1. Number of variables: 19
  2. Number of cases/rows: 76
  3. Missing data code: NA
  4. Variable list and description:
      1. plot_ID - unique identified of each plot
      2. C1_C5_L - understory canopy (<15m) in the first five sampling locations (~1m apart)
      3. C6_C10_L - understory canopy (<15m) in the second five sampling locations (~1m apart)
      4. C11_C15_L - understory canopy (<15m) in the third five sampling locations (~1m apart)
      5. C16_C20_L - understory canopy (<15m) in the fourth five sampling locations (~1m apart)
      6. C21_C25_L - understory canopy (<15m) in the fifth five sampling locations (~1m apart)
      7. C26_C30_L - understory canopy (<15m) in the sixth five sampling locations (~1m apart)
      8. C31_C35_L - understory canopy (<15m) in the seventh five sampling locations (~1m apart)
      9. C36_C40_L - understory canopy (<15m) in the eigth five sampling locations (~1m apart)
      10. C41_C45_L - understory canopy (<15m) in the ninth five sampling locations (~1m apart)
      11. C1_C5_H - overstory canopy (>20m) in the first five sampling locations (~1m apart)
      12. C6_C10_H - overstory canopy (>20m) in the second five sampling locations (~1m apart)
      13. C11_C15_H - overstory canopy (>20m) in the third five sampling locations (~1m apart)
      14. C16_C20_H - overstory canopy (>20m) in the fourth five sampling locations (~1m apart)
      15. C21_C25_H - overstory canopy (>20m) in the fifth five sampling locations (~1m apart)
      16. C26_C30_H - overstory canopy (>20m) in the sixth five sampling locations (~1m apart)
      17. C31_C35_H - overstory canopy (>20m) in the seventh five sampling locations (~1m apart)
      18. C36_C40_H - overstory canopy (>20m) in the eigth five sampling locations (~1m apart)
      19. C41_C45_H - overstory canopy (>20m) in the ninth five sampling locations (~1m apart)

## (4) DATA-SPECIFIC INFORMATION FOR: treecodes.csv

  1. Number of variables: 8
  2. Number of cases/rows: 205
  3. Missing data code: NA
  4. Variable list and description:
      1. species_ID - genus and species of identified species higher taxonomic identities (Genus, Family) or morphospecies name used to identify the species
      2. phenology_fromParkinson_fromfieldExpIfNotInParkinson - phenology status obtained by compiling information from Flora of the Andaman Islands (Parkinson, 1923), field observations and discussion with experts; E or D indicated Deciduous or Evergreen phenology where trees with over a month of leaflessness was deemed deciduous; the suffix fin was used when the phenology status was available in the flora and confirmed in the field and tmp was used when it was based on field observations alone; NA was used for higher taxonomic groups that contained species with mixed phenology within it or when the phenology of the morphospecies was not known
      3. phenologyFIN - for our study we collapsed tmp and fin suffixes and considered NA to be evergreen to reduce any Type 1 error since our null hypothesis was no change in deciduousness with logging (in the associated article)
      4. code_new - 6 letter code identical to code column in file (A)
      5. id_level_correct - the highest taxonomic level to which this morphospecies is identified - S species; G Genus; F Family; Climber liana of uknown taxonomic identity; NA no information
      6. family - Family to which the morphospecies belongs to, if known
      7. habit - tree or treelet (included); monocots include palms and cane (not included); NA includes climbers or a logging stump (not included)
      8. commercial - if the species of tree is of commercial value; Y - yes; N - no; YNC - yes but not currently logged; NA - logging stump where taxonomic identity is not known
