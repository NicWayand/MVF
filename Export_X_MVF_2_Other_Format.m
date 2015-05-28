clear all; close all; clc; smp;
tic

% This script takes Filled data in MVF (Created by
% Fill_X_Master_Format.m), formats data over a given range for a specific
% model Forcing, and saves as an user defined ascii file format, a well as
% a .mat file for reference.
%
% Nic Wayand - Jan 7th, 2013 - nicway@gmail.com
%

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
%       37=Soil Temperature A cm                    - 1(C), 2(F), 3(K)
%       38=Soil Temperature B cm                    - 1(C), 2(F), 3(K)
%       39=Soil Temperature C cm                    - 1(C), 2(F), 3(K)
%       40=Soil Temperature D cm                    - 1(C), 2(F), 3(K)
%       41=Soil Temperature E cm                    - 1(C), 2(F), 3(K)
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



%% Output Formats (Feel free to modify to define your own)
%
% output_format =
%
% 1 - FUSE Model Forcing Format
%     1    2  3  4  5     6        7              8        9          10           11     12           13
%     Year mm dd HH MM    s        PR             SW,      RDOWN,     TS,          WS,    PRES,        QS
%     1966  1  1  0  0    0.0      4.630e-05      0.000    247.200    264.250      1.000  97890.000    1.650e-03
%     1966  1  1  0 30    0.0      0.000e+00      0.000    247.300    264.288      1.010  97898.898    1.630e-03
%     1966  1  1  1  0    0.0      0.000e+00      0.000    247.300    264.308      1.018  97907.703    1.620e-03
% 
%     Where:
% 
%             PR =    precipitation rate (kg m**-2 s**-1)
%             SW =    downward shortwave radiation (W/m**2)
%          RDOWN =    downward longwave radiation (W/m**2)
%             TS =    air temperature at 2 meter height (K)
%             WS =    wind speed at 10 meter height (m/s)
%           PRES =    air pressure at 2 meter height (Pa) (NOT CORRECTED)
%             QS =    specific humidity at 2 meter height (kg/kg)

%% User Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
site = 'SNQ15';
%% Input Info
Forcing_N_in  = 'OctMar';
file_in      = strcat(site,'_Filled.mat');
%% Output Info
Forcing_N_out = 'OctMar';
file_out     = strcat(site,'_Skookum_80_08.txt');
file_out_ML  = strcat(site,'_Skookum_80_08.mat');

% Output Format
output_format = 1; 
% 1 - FUSE (Framework for Understanding Structural Errors)
% 2 - DHSVM (Distributed Hydrology Soil and Vegitation Model

% Wind speed profile scaling method (between heights)
scaling_option = 1; % 1 - Logrithmic, 2 - Exponential

%% CHECK FORMAT OPTION OUTPUT Above!!
%% CHECK WIND SPEED HEIGHTS OUTPUT IS CORRECT BELOW!!

%% UNF
% % 9697
% date_s = datenum(1996,11,25,20,0,0); % Start date
% date_e = datenum(1997,4,22,0,0,0); % End date
% UTC_offset = 0; % UTC to CEST is +1 

% % 9798
% date_s = datenum(1997,11,18,0,0,0); % Start date
% date_e = datenum(1998,4,7,22,0,0); % End date
% UTC_offset = 0; % UTC to CEST is +1

% %% CDP
% date_s = datenum(1993,8,01,1,0,0); % Start date
% date_e = datenum(2011,07,31,22,0,0); % End date
% UTC_offset = 1; % UTC to CEST is +1 

% % SNQ snotel stations 
date_s = datenum(1979,12,31,16,0,0); % Start date
date_e = datenum(2008,12,31,0,00,0); % End date
UTC_offset = 0; % UTC to CEST is +1 


% MFS
% date_s = datenum(2012,10,16,0,0,0); % Start date
% date_e = datenum(2013,05,04,0,0,0); % End date

% % SNQ 2014
% date_s = datenum(2013,10,01,0,30,0); % Start date
% date_e = datenum(2014,05,18,22,30,0); % End date

% % SNQ 2014 (Hourly)
% date_s = datenum(2013,10,01,1,0,0); % Start date
% date_e = datenum(2014,05,18,22,0,0); % End date


% SNQ 2013
% date_s = datenum(2012,10,01,0,0,0); % Start date
% date_e = datenum(2013,10,01,0,0,0); % End date

% date_s = datenum(2012,10,01,0,0,0); % Start date
% date_e = datenum(2013,06,14,09,0,0); % End date
% date_s = datenum(2013,09,09,9,0,0); % Start date
% date_e = datenum(2013,11,05,15,0,0); % End date

% % LNF
% date_s = datenum(2005,11,01,0,0,0); % Start date
% date_e = datenum(2006,01,30,21,0,0); % End date

% % % % NF97
% date_s = datenum(1995,10,01,0,0,0); % Start date
% date_e = datenum(1997,09,29,23,0,0); % End date



%% Add a Bias/offset
Adjust_Forcing = false;
Bias_1 = 1; % NO change is 1
Offset_1 = 0; % No change is 0
variable_Col_out = 9; % Column in 
if Bias_1 ~= 1 | Offset_1 ~= 0
    disp('WARNING, Bias/offset is being applied!!!!')
    pause
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Directories
homedir = 'G:\ncar\d1\data\Point_Snow_Sites\'; % REPLACE YOUR HOME DIR HERE
Indir = fullfile(homedir,'Fill',site,Forcing_N_in);
Outdir = fullfile(homedir,'Formated',site,Forcing_N_out);
if exist(Outdir,'dir') ~= 7
    mkdir(Outdir)
end
cd(Outdir)
%% Load Filled MVF data
load(strcat(Indir,'\',file_in)); 
TIME_out_all = time_builder(TIME_out);
dt_in = get_dt(TIME_out_all); %((TIME_out(2)-TIME_out(1))*24); 
sprintf('Output timestep will be %f hours',dt_in)

% Get Indices of Start and End date
I_s = find(datenum_round_off(TIME_out,'minute') == datenum_round_off(date_s,'minute')); if isempty(I_s); error('Check start date exists in Filled data'); end;
I_e = find(datenum_round_off(TIME_out,'minute') == datenum_round_off(date_e,'minute')); if isempty(I_e); error('Check end date exists in Filled data'); end;

% Create new matlab TIME
TIME_new = TIME_out_all(I_s:I_e,:);

%% Define Forcing file settings
if output_format == 1
    % FUSE MODEL Forcing Paramaters
    %   tab separated
    %   iyyy, im, id, ih, imi, dsec, PR, SW, RDOWN, T, windspd, PRES, QS,
    %   format='(i4,1x,4(i2,1x),f6.1,1x,e14.3,1x,5(f10.3,1x),e12.3)'
    MVF_out     = [2 9 11 1 5   16 4]; % MVF index
    Units_out   = [4 1 1  3 1    2 3]; % MVF units (see above options in MVF Units)
    MHeight_out = [4 4 4  4 10 4 4]; % DHSVM : Height of output varialbes (only wind speed effected)
%     MHeight_out = [4 4 4  4 7.15 4 4]; % SNQ : Height of output varialbes (only wind speed effected)
%     Time_offset_from_local = 0; % FUSE runs in Local time
    N_date_columns = 6; % Number of columns for date info
    
elseif output_format == 2
    % DHSVM MODEL FORCING INPUT
    %what: mm/dd/yyyy-HH 2m-air-temp Windspeed(50m) RH   SW    LW    P
    %unit:       -          C           m/s         %    W/m^2 W/m^2 m
    MVF_out     = [1 5  4 9 11 2];% MVF index
    Units_out   = [1 1  1 1  1 1]; % MVF units (see above options in MVF Units)
    MHeight_out = [4 50 4 4  4 4]; % m 
%     Time_offset_from_local = 0; % DSHVM runs in Local time
    N_date_columns = 4; % Number of columns for date info
else
    error('Please choose available output format')
end

%% Initialize Formated Forcing Data
N_c = (N_date_columns+length(MVF_out)); % Numbe of columns for Forcing file
N_r = (I_e-I_s+1);
Data_Forcing = NaN(N_r,N_c); 

%% Check that Required Variables exist for Forcing data and are complete, in Filled data
for cC = 1:length(MVF_out)
    N_NaNs = sum(isnan(Data_Fill(I_s:I_e,MVF_out(cC))));
    if N_NaNs > 0 % Then there are NaN's
        if N_NaNs == N_r % All data is NaN
%             alt_data_found = false; % First assume we don't have alt data source
%             
%             % Check if alternate data source exisits
%             if MVF_out(cC) == 6 % If we want wind speed at 10m
%                 N_NaNs_alt = sum(isnan(Data_Fill(I_s:I_e,5))); % Check if we have wind speed at 2m
%                 if N_NaNs_alt == 0 % Yes, we have complete wind speed at 2m
%                     alt_data_found = true;
%                     % Convert 2m to 10m Using optional profile
%                     Data_Fill(I_s:I_e,MVF_out(cC)) = scale_WS(Data_Fill(I_s:I_e,5),2,10,scaling_option);
%                     fprintf('Converting 2m wind speed to 10m wind speed using option %d \n',scaling_option);
%                     disp('Warning hard coded Xm input (instead of 2m) used here')
%                     % Keep track of previous filled data
%                     Missing_I(I_s:I_e,MVF_out(cC)) = Missing_I(I_s:I_e,5);
%                 else
%                     error('Missing enough wind speed data')
%                 end
%             elseif MVF_out(cC) == 5 % If we want wind speed at 2m
%                 N_NaNs_alt = sum(isnan(Data_Fill(I_s:I_e,6))); % Check if we have wind speed at 10m
%                 if N_NaNs_alt == 0 % Yes, we have complete wind speed at 10m
%                     lt_data_found = true;
%                     % Convert 10m to 2m Using optional profile
%                     Data_Fill(I_s:I_e,MVF_out(cC)) = scale_WS(Data_Fill(I_s:I_e,6),10,2,scaling_option);
%                     fprintf('Converting 10m wind speed to 2m wind speed using option %d \n',scaling_option);
%                     % Keep track of previous filled data
%                     Missing_I(I_s:I_e,MVF_out(cC)) = Missing_I(I_s:I_e,6);
%                 else
%                     error('Missing enough wind speed data')
%                 end
%             elseif MVF_out(cC) == 16 % If we want unadjusted pressure
%                 N_NaNs_alt = sum(isnan(Data_Fill(I_s:I_e,14))); % Check if we have an adjusted pressure
%                 if N_NaNs_alt == 0 % Yes, we have good adjusted pressure
%                     alt_data_found = true;
%                     % Convert back from adjusted to unadjusted pressure
%                     Data_Fill(I_s:I_e,MVF_out(cC)) = Press_from_sea_2_level(Data_Fill(I_s:I_e,14),site_E);
%                     fprintf('Converting sea-level presure to %3.2fm-level pressure (unadjusted) \n',site_E);
%                     % Keep track of previous filled data
%                     Missing_I(I_s:I_e,MVF_out(cC)) = Missing_I(I_s:I_e,14);
%                 else
%                     error('Missing enough wind speed data')
%                 end
%             end
%             
%             if ~alt_data_found
                fprintf('No data available for %d MVF varirable \n', MVF_out(cC))
                error('Please fill these with Fill_X_Master_Format.m')
%             end
            
            
        else % Only some NaN's exist
            fprintf('There exists %d NaN values for the %d MVF variable \n',N_NaNs,MVF_out(cC))
            error('Please fill these with Fill_X_Master_Format.m')
        end
    else % There are no NaN's present, data is complete
        % Do nothing
    end
    clear N_NaNs
end

%% Check that output variables are at the correct height

for cC = 1:length(MVF_out)
    if MHeight_MVF(MVF_out(cC)) ~= MHeight_out(cC); 
        if MVF_out(cC) == 5 || MVF_out(cC) == 6 % Need to adjust windspeed height (based on method above)
            Data_Fill(I_s:I_e,MVF_out(cC)) = scale_WS(Data_Fill(I_s:I_e,MVF_out(cC)),MHeight_MVF(MVF_out(cC)),MHeight_out(cC),scaling_option);
            fprintf('Converting %d m wind speed to %d m wind speed using option %d \n',MHeight_MVF(MVF_out(cC)),MHeight_out(cC),scaling_option);
        end
    end  
end

%% Convert Time to output format
if output_format == 1
    % FUSE MODEL Forcing Paramaters
    TIME_out_all = time_shift(TIME_out_all,UTC_offset);
    Data_Forcing(:,1:5) = TIME_out_all(I_s:I_e,1:5);
    Data_Forcing(:,6) = 0; % Seconds hardcoded to zero - NIC
elseif output_format == 2
    % DHSVM [2 3 1 4]
    TIME_out_all = time_shift(TIME_out_all,UTC_offset);
    Data_Forcing(:,1:4) = [TIME_out_all(I_s:I_e,2) TIME_out_all(I_s:I_e,3) TIME_out_all(I_s:I_e,1) TIME_out_all(I_s:I_e,4)];
end

%% Grab time series of Forcing Variables and Convert Units if required
UNIT_Forcing = ones(1,length(Units_out)); % Assumed all units have been convered to MVF format in Data_Fill
for cC = 1:length(MVF_out)
    
    if Units_out(cC) == 1 % Forcing units are MVF units, no conversion needed
        Data_Forcing(:,N_date_columns+cC)  = Data_Fill(I_s:I_e,MVF_out(cC));
    else % Some conversion is needed
        % What type of variable is it?
        if sum(MVF_out(cC) == [1 22 31 42 43]) % Air or surface temperature
            if Units_out(cC) == 2     % Convert C to F
                Data_Forcing(:,N_date_columns+cC)  = C_to_F(Data_Fill(I_s:I_e,MVF_out(cC))); UNIT_Forcing(cC) = 2;
            elseif Units_out(cC) == 3 % Convert C to K
                Data_Forcing(:,N_date_columns+cC)  = C_to_K(Data_Fill(I_s:I_e,MVF_out(cC))); UNIT_Forcing(cC) = 3;
            else
                error('Unknown unit specified, check list')
            end  
        elseif sum(MVF_out(cC) == [2 3 20 23 24 44]) % Precip, Snow Depth
            if Units_out(cC) == 2     % Convert m to mm
                Data_Forcing(:,N_date_columns+cC)  = m_to_mm(Data_Fill(I_s:I_e,MVF_out(cC))); UNIT_Forcing(cC) = 2;
            elseif Units_out(cC) == 3 % Convert m to in
                Data_Forcing(:,N_date_columns+cC)  = m_to_in(Data_Fill(I_s:I_e,MVF_out(cC))); UNIT_Forcing(cC) = 3;
            elseif Units_out(cC) == 4 % Convert m to kg m**-2 s**-1
                Data_Forcing(:,N_date_columns+cC)  = m_to_kg_per_m2_s(Data_Fill(I_s:I_e,MVF_out(cC)),dt_in); UNIT_Forcing(cC) = 4;
            else
                error('Unknown unit specified, check list')
            end
        elseif sum(MVF_out(cC) == [4 32]) % Air Moisture content
            if Units_out(cC) == 2     % Convert RH to fraction
                Data_Forcing(:,N_date_columns+cC)  = RH_to_fraction(Data_Fill(I_s:I_e,MVF_out(cC))); UNIT_Forcing(cC) = 2;
            elseif Units_out(cC) == 3     % Convert RH to Mixing Ratio Q (kg/kg)
                Data_Forcing(:,N_date_columns+cC)  = RH_to_Q(Data_Fill(I_s:I_e,MVF_out(cC)),...
                    Data_Fill(I_s:I_e,1), Data_Fill(I_s:I_e,(16)) ); UNIT_Forcing(cC) = 3;
            else
                error('Unknown unit specified, check list')
            end   
        elseif sum(MVF_out(cC) == [14 16]) % Air Pressure
            if Units_out(cC) == 2     % Convert hPa to Pa
                Data_Forcing(:,N_date_columns+cC)  = hPa_to_Pa(Data_Fill(I_s:I_e,MVF_out(cC))); UNIT_Forcing(cC) = 2;
            elseif Units_out(cC) == 3     % Convert hPa/mb to kPa 
                Data_Forcing(:,N_date_columns+cC)  = hPa_to_kPa(Data_Fill(I_s:I_e,MVF_out(cC))); UNIT_Forcing(cC) = 3;
            else
                error('Unknown unit specified, check list')
            end 
           
        else
            error('need other options')
        end


    end
    
    
end
fprintf('All units converted to format %d \n', output_format)

% Adjust for bias or offset in code
if Adjust_Forcing
    Data_Forcing(:,variable_Col_out) = Data_Forcing(:,variable_Col_out) .* Bias_1 + Offset_1;
    fprintf('Adjusted variable %d by an bias of %f and an offset of %f\n',variable_Col_out,Bias_1,Offset_1)
end

% Checks
if sum(sum(isnan(Data_Forcing))) > 0
    error('Warning, NaNs still remain in the forcing file')
end

% Export to ascii format as specified
if output_format == 1 
    % FUSE
    % format='(i4,1x,4(i2,1x),f6.1,1x,e14.3,1x,5(f10.3,1x),e12.3)'
    % 1966 10 15 23 30    0.0      0.000e+00      0.000    323.900    276.299      1.382 100133.398    4.860e-03
    fid = fopen(fullfile(Outdir,file_out),'w');
    fprintf(fid,'%4i %2i %2i %2i %2i %6.1f %14.3e %10.3f %10.3f %10.3f %10.3f %10.3f %12.3e \n',Data_Forcing(:,1:13)');
    fclose(fid);

elseif output_format == 2
    % DHSVM
    fid = fopen(fullfile(Outdir,file_out),'w');
    for crow = 1:length(Data_Forcing)
        %what: mm/dd/yyyy-HH 2m-air-temp Windspeed(50m) RH   SW    LW    P
        fprintf(fid,'%02i/%02i/%02i-%02i ',Data_Forcing(crow,1:4)'); % write time
        if crow ~=length(Data_Forcing)
            fprintf(fid,'%f %f %f %f %f %.10f\n',Data_Forcing(crow,5:10)'); % write data
        elseif crow ==length(Data_Forcing)
            fprintf(fid,'%f %f %f %f %f %.10f',Data_Forcing(crow,5:10)'); % write data
        end
    end

%     fprintf(fid,'%2i/%02i/%02i-%02i %6.3f %6.3f %3.3f %6.3f %6.3f %4.10f \n',Data_Forcing(:,1:10)');
    fclose(fid);
else
    error('not ready')
end

% Save Matlab version of Forcing data
Missing_I_Forcing = Missing_I(I_s:I_e,MVF_out);
Missing_I_Forcing = [zeros(size(Missing_I_Forcing,1),6) Missing_I_Forcing]; % add zeros to date columns
eval(['save ' file_out_ML ' Data_Forcing TIME_new site_E Missing_I_Forcing site_E'])

disp(' ')
fprintf('Finished Exporting data for forcing %s \n',Forcing_N_out)
