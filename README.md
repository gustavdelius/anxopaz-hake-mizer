# hake-mizer

Repository to work the hake model with Mizer for the FRESCO project. The scripts must be run in the order in which they are indicated for a correct operation:  

1-Bio_Pars.R: Allows to create the model with the biological data obtained for the southern hake stock in ICES from a series of files available in the folder "scripts" where each biological process is modeled or read directly from SS3. These scripts contain the input and output information for each process and where they are available. In addition, results plots are generated.  

2-Catch.R: Obtaining the catch information from the data available in the “data” folder transformed in such a way that they can be used later in the adjustment.  

3-MIZER.R: The model for hake is built from the biological information, catches and data obtained from the current assessment for this stock and the selectivity is adjusted from a C++ function contained in the "TMB" folder.  

