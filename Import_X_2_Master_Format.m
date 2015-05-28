clear all; close all; clc;
tic

%% Creating Forcing data using the Master Variable Format 
% Nic Wayand - May 27th, 2015 - nicway@gmail.com
%
% This script is used to import various formats of raw data into the Master
% Variable Format (MVF), which allows for the following Quality Control, Filling, and
% Formating scripts to use the same format of data between different sites
% or years. This is the first of 4 scripts used to properly take met data
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
%		42=Air Temperature C						- 1(C), 2(F), 3(K)
%		43=Air Temperature D						- 1(C), 2(F), 3(K)
%		44=Height of CSAT (from ground)				- 1(m), 2(mm), 3(in), 4(kg m**-2 s**-1), 5(cm), 6(k-ft)
%		45=CSAT Azimuth from true north		    	- 1(degree)
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

MVF_N_vars = 70;  % Number of total variable columns to allow. This is flexible.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% User Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Primary data 
% Snoqualmie 2015 data  (Oct-March) 5/27/2015
site           = 'SNQ15'; % Site folder name
Forcing_N      = 'OctMar'; % Forcing project
site_E         = 917; % m
appending_data = false; % false - creates new mat file, true - adds to existing file (throws error if asked to overwrite any previous column with data)

% SNQ data (UW and NWAC tower) 
file_in        = ('halfhour_SPR2015.mat'); header_lines = NaN; date_format = 3; % 1- SNQ, 2- MFS, 3- Matlab format
TIME_IN = half_15.timePST(:,7); % Define name of matlab time variable Nx7
DATA_IN = half_15.data; % Define name of matlab data variable
NA             = -99; % Unused variable (use in MVF variable if you don't want to import said column)
time_offset    = 0; % Hours offset from input time to local standard time

% Index of MVF refers to columns in the DATA_IN, while the value within the variable MVF
% refers to the index as defined above in the Master Variable Formate.
% raw data index                           NWAC ---------------------------------- NWAC  UW ------------------------------------------------
               %    1   2   3   4   5   6   7   8   9   10   11      12   13   14   15   16   17   18   19   20   21   22   23   24   25   26   27   28   29   30   31   32   33   34   35   36   37   38   39   40   41   42   43   44   45   46   47   48   49  50   51   52   53   54   55   56   57   58   59    60    61    62   63   64   65   66   67   68   69   70   71   72   73   74    75    76    77     78    79    80    81    82   83   84   85   86    87   88   89   90    % Indice from input data 
MVF            = [-99 -99 -99 -99 -99 -99  31  32   2   14    5      15  -99  -99   20   23  -99  -99  -99  -99  -99  -99  -99  -99  -99  -99  -99  -99  -99  -99   46  -99  -99  -99  -99  -99  -99   42  -99  -99  -99  -99  -99    9   10  -99  -99  -99  -99  11   12  -99  -99   19   13  -99   22  -99  -99    33    34    35  -99  -99  -99    1  -99  -99  -99    4  -99   43   44   45    38    39    40     37    56    57    58    59  -99  -99  -99  -99    33   51  -99  -99 ]; % Indices of MVF (see above) within the input file          
%                               
UNIT           = [  1   1   1   1   1   1   1   1   2    1    1       1    1    1    1    2    1    1    1    1    1    1    1    1    1    1    1    1    1    1    3    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1   1    1    1    1    1    1    1    1    1    1     1     1     1    1    1    1    1    1    1    1    2    1    1    1    1     1     1     1      1     1     1     1     1    1    1    1    1     1    1    1    1 ]; % Units of the input variable as defined above
MHeight        = [-99 -99 -99 -99 -99 -99   4   4   7.15 6    7.15    6  -99  -99    0    0  -99  -99  -99  -99  -99  -99  -99  -99  -99  -99  -99  -99  -99  -99    2  -99  -99  -99  -99  -99  -99    6  -99  -99  -99  -99  -99    6    6  -99  -99  -99  -99   6    6  -99  -99    6    6  -99    0  -99  -99  -.06  -.06  -.06  -99  -99  -99    6  -99  -99  -99    6  -99    5    0    0 -0.10 -0.19 -0.405 -0.02 -0.02 -0.10 -0.19 -0.36  -99  -99  -99  -99 -0.08    0  -99  -99 ]; % Height of measurements in m, corresponding to A, B, C in MVF

% 
% % If appending data, comment out above and uncomment block below
% %%%%%%%%%%%%% Appending data (snowpit.mat)
% site           = 'SNQ15'; % Site folder name
% Forcing_N      = 'OctMar'; % Forcing project
% site_E         = 917; % m
% appending_data = true; % false - creates new mat file, true - adds to existing file (throws error if asked to overwrite any previous column with data)
% 
% % % SNQ data (24 6am observations) 
% file_in        = ('Snoqualmie_Daily_2014-2015.mat'); header_lines = NaN; date_format = 3; % 1- SNQ, 2- MFS, 3- Matlab format
% TIME_IN = snowpit.timePST; % Define name of matlab time variable Nx7
% DATA_IN = snowpit.data; % Define name of matlab data variable
% NA             = -99; % Unused variable (use in MVF variable if you don't want to import said column)
% time_offset    = 0; % Hours offset from input time to local time
% % %                 1  2  3  4  5  6    7  8   9   10  11
% MVF            = [-99 48 24 49 50 51   52 53  54  -99  55];
% UNIT           = [  1  5  5  5  1  1    1  6   5    1   3];
% MHeight        = [-99  0  0  0  0  0 -0.2  0   0  -99   0];
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Directories and output files
homedir = 'G:\ncar\d1\data\Point_Snow_Sites\'; % REPLACE YOUR HOME DIR HERE
rawdir = fullfile(homedir,'Raw',site,Forcing_N);

file_out = strcat(site,'_RAW.mat');
file_save = fullfile(rawdir,file_out); 
if exist(rawdir,'dir') ~= 7
    mkdir(rawdir)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If we are appending a previous mat file
if appending_data
    load(file_save)
end

%% Checks

if size(unique(MVF(MVF~=NA))) ~= size(MVF(MVF~=NA))
    error('Check MVF, there are repeated variables')
end

% Set date format of input data
if date_format == 1
    % Input date format
    % {'YEAR','DOY','HOUR'}
    TIME_columns = 3;
    %% Column numbers for YYYY, DOY, and HH (hour)
    YYYY_in = 1;
    DOY_in = 2; 
    HH_in = 3; HH_adj = 1./100; % to HH format
    MM_val_in = 00;
    SS_val_in = 00;
elseif date_format == 2
    % Input date format
    % {'2001-09-27 20:24:00';}
    TIME_columns = 1;
%     TIME_in = datenum(time_in); % too slow here
elseif date_format == 3
    % Already in matlab format. Do nothing

elseif date_format == 4 % CDEC copy and past 
    % 06/03/1996 00:00        65
    TIME_columns = 1;
    
else
    error('Please select avilable option')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   DO WORK  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Import ascii file(s)
cd(rawdir)

% Input Raw data
if date_format == 1 % You will have to update below to suite your import code format
    data_all = dlmread(file_in,',',header_lines, 0); 
    data_in = data_all(:,TIME_columns+1:end); 
    time_in = data_all(:,1:TIME_columns); clear data_all
    
    % Get first and last date of input data, and input time step (assumed constant) 
    [cmonth, cday] = jul2greg(time_in(1,YYYY_in),time_in(1,DOY_in)); % Get Month and Day
    Time_s = datenum(time_in(1,YYYY_in),cmonth,cday,time_in(1,HH_in)*HH_adj,MM_val_in,SS_val_in); % build time
    [cmonth, cday] = jul2greg(time_in(2,YYYY_in),time_in(2,DOY_in)); % Get Month and Day
    Time_2 = datenum(time_in(2,YYYY_in),cmonth,cday,time_in(2,HH_in)*HH_adj,MM_val_in,SS_val_in); % build time
    [cmonth, cday] = jul2greg(time_in(end,YYYY_in),time_in(end,DOY_in)); % Get Month and Day
    Time_e = datenum(time_in(end,YYYY_in),cmonth,cday,time_in(end,HH_in)*HH_adj,MM_val_in,SS_val_in); % build time
    dt_in = floor((Time_2-Time_s)*24); % in hours
    
    if ~appending_data 
        % Make new matlab time series for MVF
        TIME_out_all = time_builder(Time_s,Time_e,dt_in);
        TIME_out = TIME_out_all(:,7); 
    else
        TIME_out_all = time_builder(TIME_out(1),TIME_out(end),dt_in);
    end
    
elseif date_format == 2 % You will have to update below to suite your import code format
    data_in = dlmread(file_in,',',header_lines, TIME_columns); 
    data_tree = importdata(file_in);
    time_in = data_tree.textdata(header_lines+1:end,TIME_columns); clear data_tree
    % Need to grab text date in a better way (using importdata for now)
    % time_in = csvread(file_in, header_lines, 0, [header_lines 0 length(data_in) 0]);
    % format_time = '%s %s'; 
    % time_in = textscan(file_in,format_time

    Time_s = datenum(time_in(1,1));
    Time_e = datenum(time_in(end,1));
    dt_in = (datenum(time_in(2,1)) - datenum(time_in(1,1)))*24; % in hours 
    
    if ~appending_data 
        % Make new matlab time series for MVF
        TIME_out_all = time_builder(Time_s,Time_e,dt_in);
        TIME_out = TIME_out_all(:,7); 
    else
        TIME_out_all = time_builder(TIME_out(1),TIME_out(end),dt_in);
    end
    
elseif date_format == 3
    load(file_in)
    
    NewTimeALL = datenum_round_off(TIME_IN,'minute'); % Round off matlab date errors
    time_in = NewTimeALL(:,7); % Grab serial dates
    data_in = DATA_IN; % Get data
   
    if ~appending_data % TIME out = TIME in
        TIME_out     = time_in;
        TIME_out_all = time_builder(TIME_out);
    else % Use original TIME out for TIME out
        % do nothing, orignal TIME_out an TIME_out_all were loaded above
    end
    
elseif date_format == 4 % You will have to update below to suite your import code format
    data_tree = importdata(file_in,',',0);
    time_in = data_tree.textdata(header_lines+1:end,TIME_columns); 
    data_in = double(data_tree.data);
    
    Time_s = datenum_round_off(datenum(time_in(1,1)),'minute');
    Time_e = datenum_round_off(datenum(time_in(end,1)),'minute');
    disp('Hard coded timestep to 1 hour!!!!!!!!!!!!!!!!!!!!!!!!!')
    dt_in = 1; %floor(datenum_round_off(datenum(time_in(2,1)),'minute') - datenum_round_off(datenum(time_in(1,1)),'minute'))*24.00000); % in hours 
    
    if ~appending_data 
        % Make new matlab time series for MVF
        TIME_out_all = time_builder(Time_s,Time_e,dt_in);
        TIME_out = TIME_out_all(:,7); 
    else
        TIME_out_all = time_builder(TIME_out(1),TIME_out(end),dt_in);
    end
end
%% CHECKS

% Check data_in has the same columns as specified in MVF
if size(data_in,2) ~= length(MVF)
    error('Columns in MVF must match input data columns (data_in)')
end

% Initialize
Nts = size(data_in,1);
if ~appending_data
    Data_Raw = nan(size(TIME_out,1), MVF_N_vars); 
    MHeight_MVF = nan(1,MVF_N_vars);
    
    % Build MVF index of Measurment heights (from MHeight)
    MHeight_MVF(MVF(MVF~=NA)) = MHeight(MVF~=NA);
end

% Round time to nearest minute to allow different data sets to be combined
TIME_out = datenum_round_off(TIME_out,'minute');
disp('')
disp('Rounding TIME to nearest minute')
disp('')

%% Add Raw data into Master Variable Format (MVF)
disp('')
disp('Add Raw data into Master Variable Format (MVF)')
disp('')
% Must check each time step (may be duplicates or missing in input)
Unfounddates = 0; 
for cts = 1:Nts 
    % Make Matlab date for each input time step
    if date_format == 1
        if time_in(cts,DOY_in) >= 1 && time_in(cts,DOY_in) <=366 % Check good DOYTIME_out
            [cmonth, cday] = jul2greg(time_in(cts,YYYY_in),time_in(cts,DOY_in));
            cTime = datenum(time_in(cts,YYYY_in),cmonth,cday,time_in(cts,HH_in)*HH_adj,MM_val_in,SS_val_in);
        else
            cTime = -9999; % Our value for a bad date
        end
    elseif date_format == 2
        cTime = datenum(time_in(cts));
    elseif date_format == 3
        cTime = time_in(cts);
    elseif date_format == 4
        try
            cTime = datenum(time_in(cts));
        catch
            fprintf('Found bad time --> %s, which we throw away\n',char(time_in(cts)))
            cTime = -99999; % Set negative to skip for loop below
        end
    end
    
    % Find Index in output TIME_out
    if cTime > 0 % Positive if a good date was provided
        I_out = find(TIME_out == cTime); 

        if ~isempty(I_out) % Atleast one date was found

            if numel(I_out) > 1 % Many were found
                sprintf('Found %i records for date %s, assuming first record is correct',numel(I_out),datestr(cTime))
                I_out = I_out(1); 
            end

            % CURRENTLY SLOW, REDO FOR SPEED (need to handel unused columns)
            %Loop through each column in data_in
            %If appending, only replace values for which Data_Raw is nan
            %and datain has non-Nan values.
            for cC = 1:size(data_in,2) 
                
               %   % import it %        % Current value is a NaN %     % New values is non-NaN %
               if MVF(cC) ~= -99  %| ( isnan(Data_Raw(I_out,MVF(cC))) & ~isnan(data_in(cts,cC)) )
                  Data_Raw(I_out,MVF(cC)) = data_in(cts,cC); 
               end
               
            end
        else
            Unfounddates(cts) = cTime;
        end
        clear I_out cTime
        
    else
        % It was a bad date, do nothing
    end
end
disp('Done')

% Check fo runfound dates

if sum(Unfounddates) > 0
    disp('WARNING!!! check Unfounddates to see what timesteps was not imported!')
end

%% Format to MVF Units (1)
disp('')
disp('Checking and Converting Units')
disp('')
%Loop through each column in MFV
UNIT_new = UNIT; % To track changes
for cC = 1:size(data_in,2); 
    if UNIT(cC) ~=1 % Not in MVF unit (1), then convert
        % What type of variable is it?
        if sum(MVF(cC) == [1 22 31 37 38 39 40 41 42 43 51 52 ]) % Air or surface temperature
            if UNIT(cC) == 2     % Convert F to C
                Data_Raw(:,MVF(cC)) = F_to_C(Data_Raw(:,MVF(cC))); UNIT_new(cC) = 1;
            elseif UNIT(cC) == 3 % Convert K to C
                Data_Raw(:,MVF(cC)) = K_to_C(Data_Raw(:,MVF(cC))); UNIT_new(cC) = 1;
            else
                error('Unknown unit specified, check list')
            end  
        elseif sum(MVF(cC) == [2 3 20 23 24 44 47 48 49 53 54 55 63]) % Precip, Snow Depth, lysimeter, SWE
            if UNIT(cC) == 2     % Convert mm to m
                Data_Raw(:,MVF(cC)) = mm_to_m(Data_Raw(:,MVF(cC))); UNIT_new(cC) = 1;
            elseif UNIT(cC) == 3 % Convert in to m
                Data_Raw(:,MVF(cC)) = in_to_m(Data_Raw(:,MVF(cC))); UNIT_new(cC) = 1;
            elseif UNIT(cC) == 5 % Convert cm to m
                Data_Raw(:,MVF(cC)) = cm_to_m(Data_Raw(:,MVF(cC))); UNIT_new(cC) = 1;
            elseif UNIT(cC) == 6 % Convert 1000's ft to m
                Data_Raw(:,MVF(cC)) = kft_to_m(Data_Raw(:,MVF(cC))); UNIT_new(cC) = 1;
            else
                error('Unknown unit specified, check list')
            end
        elseif sum(MVF(cC) == [4 32]) % Air Moisture content
            if UNIT(cC) == 2     % Convert fraction to RH
                Data_Raw(:,MVF(cC))  = fraction_to_RH(Data_Raw(:,MVF(cC))); UNIT_new(cC) = 1;
            elseif UNIT(cC) == 3     % Convert Mixing Ratio Q (kg/kg) to RH
                Data_Raw(:,MVF(cC))  = Q_to_RH(Data_Raw(:,MVF(cC)),...
                    Data_Raw(:,1), Data_Raw(:,(14)) ); UNIT_new(cC) = 1;
            else
                error('Unknown unit specified, check list')
            end   
        elseif sum(MVF(cC) == [14 16 46]) % Air Pressure
            if UNIT(cC) == 2     % Convert Pa to hPa/mb
                Data_Raw(:,MVF(cC))  = Pa_to_hPa(Data_Raw(:,MVF(cC))); UNIT_new(cC) = 1;
            elseif UNIT(cC) == 3     % Convert kPa to hPa/mb
                Data_Raw(:,MVF(cC))  = kPa_to_hpa(Data_Raw(:,MVF(cC))); UNIT_new(cC) = 1;
            else
                error('Unknown unit specified, check list')
            end 
        elseif sum(MVF(cC) == [5 6]) % Wind Speed
            if UNIT(cC) == 2     % Convert mph to m/s
                Data_Raw(:,MVF(cC))  = mph_to_mps(Data_Raw(:,MVF(cC))); UNIT_new(cC) = 1;
            else
                error('Unknown unit specified, check list')
            end
        elseif sum(MVF(cC) == [56 57 58 59 60 61]) % Soil Moisture
            if UNIT(cC) == 2     % ???
                error('need to make conversion scripts')
%                 Data_Raw(:,MVF(cC))  = mph_to_mps(Data_Raw(:,MVF(cC))); UNIT_new(cC) = 1;
            else
                error('Unknown unit specified, check list')
            end
        else
            error('need other options')
        end

    end
   
    
end

% Check all Units are now in MVF Units (i.e. 1)
if any(UNIT_new > 1)
    error('Some of UNIT_new is still not converted to MVF units!')
end
disp('Done')

%% Change time zone to Local Standard Time
if time_offset > 0
    disp('Shifting Time to Local Standard Time')
    TIME_out_all = time_shift(TIME_out_all,time_offset);
    TIME_out = TIME_out_all(:,7);
    disp('Done')
end

%% Save output
if ~appending_data
    eval(['save ' file_save ' Data_Raw TIME_out MHeight_MVF site_E'])
else    
    disp('Please check Data_Raw output to see if data as appended correctly')
    disp(' ')
    pause
    disp('Enter 1 to save this appended data, Or 2 to revert to previous data')
    ui1 = user_inp(cellstr('...'), 1, 2);
    if ui1 == 1
        eval(['save ' file_save ' Data_Raw TIME_out MHeight_MVF site_E'])
        fprintf('New data saved as %s \n',file_save)
    elseif ui1 == 2
        % do nothing, file_save remains as previous data
    end

end
    
%% END
toc