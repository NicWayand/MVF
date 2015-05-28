# MVF
Matlab scripts to import, quality control, fill and export formated forcing data for point snow model simualtions.

% These scripts are used to import various formats of raw data into the Master
% Variable Format (MVF), which allows for the following Quality Control, Filling, and
% Formating scripts to use the same format of data between different sites
% or years. There are 4 scripts used to properly take met data
% from a raw state to a form used to drive any given model.
%
% Order of script to Run.
%
% 1) Import_X_2_Master_Formate.m         - Imports raw data to MVF
% 2) QC_X_Master_Formate.m               - Interactivily QC data
% 3) Fill_X_Master_Format.m              - Fill or replace missing data
% 4) Export_X_MVF_2_Other_Format.m       - Format select variables for input
%
% Search for "REPLACE YOUR HOME DIR HERE" to update paths
%
% Suggested file structure:
%
% ~/home/data/
%            /raw/
%                /site1/
%                      /Forcingfile_2000/
%                      /Forcingfile_All_years/
%                /site2/
%            /QC/
%                /site1/
%                      /Forcingfile_2000/
%                      /Forcingfile_All_years/
%                /site2/
%            /Fill/
%                /site1/
%                      /Forcingfile_2000/
%                      /Forcingfile_All_years/
%                /site2/
%            /Formated/
%                /site1/
%                      /Forcingfile_2000/
%                      /Forcingfile_All_years/
%                /site2/
%
% This structure allows the user to easily keep track of the history of a
% given forcing set. In addition, an individual variable may be updated by
% it self (i.e. additional filling from a new source) without remaking all
% other variables.


% The MVF is described below. 
%%       Master Variable Format (MVF)
%
%       # = Variable                   - Unit #
%       A - Primary, B - Secondary (For multiple observations for a given variable)
%       
%       1=Air temperature A                         - 1(C), 2(F), 3(K) 
%       2=Incremental Precipitation A               - 1(m), 2(mm), 3(in), 4(kg m**-2 s**-1), 5(cm), 6(k-ft)
%       3=Accumulated Precipitaiton A               - 1(m), 2(mm), 3(in), 4(kg m**-2 s**-1), 5(cm), 6(k-ft)
%       4=Air Moisture Content A                    - 1(RH %), 2(RH Fraction), 3(Q kg/kg)
%       5=Scalar Wind Speed A                       - 1(m/s), 2(mph)
%       6=Scalar Wind Speed B                       - 1(m/s), 2(mph)
%       7=Wind Direction at A                       - 1(degrees)
%       8=Wind Direction at B                       - 1(degrees)
%       9=Downward Solar Radiation                  - 1(W/m^2), 2(MJ/timestep)
%       10=Upward Soloar Radiation                  - 1(W/m^2), 2(MJ/timestep) 
%       11=Downward Terrestrial Rad.                - 1(W/m^2), 2(MJ/timestep)
%       12=Upward Terrestrial Rad.                  - 1(W/m^2), 2(MJ/timestep)
%       13=Net Radiation                            - 1(W/m^2), 2(MJ/timestep)
%       14=Total Pressure adjusted to sea-level		- 1(mb/hPa), 2(Pa), 3(kPa)
%       15=Rads (Sunny measure)                     - 1(1-10, unitless)
%       16=Total Pressure unadjusted A              - 1(mb/hPa), 2(Pa), 3(kPa)
%       17=Dew Point                                - 1(C)
%       18=Wetness sensor                           - 1(mV), 2(V)
%       19=Alebdo                                   - 1(fraction)
%       20=Snow Depth A                             - 1(m), 2(mm), 3(in), 4(kg m**-2 s**-1), 5(cm), 6(k-ft)
%       21=Snow Depth QC value                      - 1(user defined)
%       22=Surface Temperature A                    - 1(C), 2(F), 3(K)
%       23=Lysimeter outflow                        - 1(m), 2(mm), 3(in), 4(kg m**-2 s**-1), 5(cm), 6(k-ft)
%       24=Snow Depth B (6am snow-stake manual)     - 1(m), 2(mm), 3(in), 4(kg m**-2 s**-1), 5(cm), 6(k-ft)
%       25=Sensible Heat Flux                       - 1(w/m^2), 2(MJ/timestep)
%       26=Latent Heat Flux                         - 1(w/m^2), 2(MJ/timestep)
%       27=Soil Reflectometer Output Period 10cm    - 1(usec) 
%       28=Soil Wetness Fraction (x100) 10cm        - 1(%)
%       29=Soil Reflectometer Output Period  15cm   - 1(usec)
%       30=Soil Wetness Fraction (x100)  15cm       - 1(%)
%       31=Air temperature B                        - 1(C), 2(F), 3(K)
%       32=Air Moisture Content B                   - 1(RH %), 2(RH Fraction), 3(Q kg/kg)
%       33=Soil Heat Flux  (1)                      - 1(w/m^2), 2(MJ/timestep)
%       34=Soil Heat Flux  (2)                      - 1(w/m^2), 2(MJ/timestep)
%       35=Soil Heat Flux  (3)                      - 1(w/m^2), 2(MJ/timestep)
%       36=Battery Voltage                          - 1(volts)
%       37=Soil Temperature A                       - 1(C), 2(F), 3(K)
%       38=Soil Temperature B                       - 1(C), 2(F), 3(K)
%       39=Soil Temperature C                       - 1(C), 2(F), 3(K)
%       40=Soil Temperature D                       - 1(C), 2(F), 3(K)
%       41=Soil Temperature E                       - 1(C), 2(F), 3(K)
%		    42=Air Temperature C					            	- 1(C), 2(F), 3(K)
%		    43=Air Temperature D					            	- 1(C), 2(F), 3(K)
%		    44=Height of CSAT (from ground)		      		- 1(m), 2(mm), 3(in), 4(kg m**-2 s**-1), 5(cm), 6(k-ft)
%		    45=CSAT Azimuth from true north		        	- 1(degree)
%       46=Total Pressure unadjusted B              - 1(mb/hPa), 2(Pa), 3(kPa)
%       47=Snow Water Equivelent A                  - 1(m), 2(mm), 3(in), 4(kg m**-2 s**-1), 5(cm), 6(k-ft)
%       48=24hr Snow Depth change 6am               - 1(m), 2(mm), 3(in), 4(kg m**-2 s**-1), 5(cm), 6(k-ft)
%       49=24hr SWE change 6am                      - 1(m), 2(mm), 3(in), 4(kg m**-2 s**-1), 5(cm), 6(k-ft)
%       50=24hr Density 6am                         - 1(%), 2(fraction)
%       51=Surface Temperature B                    - 1(C), 2(F), 3(K)
%       52=Snow Layer Temperature A                 - 1(C), 2(F), 3(K)
%       53=Freazing level A                         - 1(m), 2(mm), 3(in), 4(kg m**-2 s**-1), 5(cm), 6(k-ft)
%       54=RamPen Test A                            - 1(m), 2(mm), 3(in), 4(kg m**-2 s**-1), 5(cm), 6(k-ft)
%       55=Incremental Precipitation B              - 1(m), 2(mm), 3(in), 4(kg m**-2 s**-1), 5(cm), 6(k-ft)
%       56=Soil Moisture A                          - 1(VWC,-)
%       57=Soil Moisture B                          - 1(VWC,-)
%       58=Soil Moisture C                          - 1(VWC,-)
%       59=Soil Moisture D                          - 1(VWC,-)
%       60=Soil Moisture E                          - 1(VWC,-)
%       61=Soil Moisture F                          - 1(VWC,-)
%       62=Surface Temperature C                    - 1(C), 2(F), 3(K)
%       63=Lysimeter outflow B                      - 1(m), 2(mm), 3(in), 4(kg m**-2 s**-1), 5(cm), 6(k-ft)

