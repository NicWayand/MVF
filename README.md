# MVF (Master Variable Format)
Matlab scripts to import, quality control, fill and export formated forcing data for point snow model simualtions.

 These scripts re used to import various formats of raw data into the Master
 Variable Format (MVF), which allows for the following Quality Control, Filling, and
 Formatting scripts to use the same format of data between different sites
 or years. There are 4 scripts used to properly take met data
 from a raw state to a form used to drive any given model.

 Order of script to Run.

 1) Import_X_2_Master_Formate.m         - Imports raw data to MVF
 2) QC_X_Master_Formate.m               - Interactivily QC data
 3) Fill_X_Master_Format.m              - Fill or replace missing data
 4) Export_X_MVF_2_Other_Format.m       - Format select variables for input

 Search for "REPLACE YOUR HOME DIR HERE" to update paths

 Suggested file structure:

 ~/home/data/
            /raw/
                /site1/
                      /Forcingfile_2000/
                      /Forcingfile_All_years/
                /site2/
            /QC/
                /site1/
                      /Forcingfile_2000/
                      /Forcingfile_All_years/
                /site2/
            /Fill/
                /site1/
                      /Forcingfile_2000/
                      /Forcingfile_All_years/
                /site2/
            /Formated/
                /site1/
                      /Forcingfile_2000/
                      /Forcingfile_All_years/
                /site2/

 This structure allows the user to easily keep track of the history of a
 given forcing set. In addition, an individual variable may be updated by
 it self (i.e. additional filling from a new source) without remaking all
 other variables.
