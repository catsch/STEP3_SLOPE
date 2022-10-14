# STEP3_SLOPE

In order to perform DM for CHLA you should
1. STEP 1 : Estimate the DARK OFFSET of the CHLA

- https://github.com/catsch/dark_offset_chla on one core
  or
- https://github.com/qjutard/dark_offset_chla on several cores

2. STEP 2 : Correct the profile from the DARK and estimate the quenching correction (Xing 2018, Terrats 2020, Sackmann 2008)

- https://github.com/catsch/STEP2_QUENCHING

3. STEP3 : Estimate the F490 factor (Slope factor that converts your Fluorescence profile corrected of the DARK and the Quenching into a CHLA profile) 

  => A working directory containing Argo B files and C files (OUTPUTS of the STEP2)

  => coriolis_CHLA.list is the list of all the WMO floats for which you want to estimate your slope 

  -lance_STEP3 is my launching script :

      -It defines the pathway to the working directory

      -It reads the coriolis_CHLA.list to know the WMO

      -It reads the working directory to build the list of the Bfiles : "liste_all_B"

  -SLOPE_ESTIMATION.R is preparing the entries for Xing2011.R 

  -Xing2011.R is estimating the F490/SLOPE  https://doi.org/10.17882/86384)

  -RunningFilter.R is used to filter (mean, median) the data


  The outputs are files containing the SLOPE estimation for the profiles and QC (1,2 is good) 

