clear all; close all; clc; smp;
tic

% This script takes Raw data in MVF (Created by
% Import_X_2_Master_Format.m), performs Quality Control on each variable,
% and saves a X_QC.mat, Wher X is the 3 letter station code
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
%       62=Surface Temperature C                    - 1(C), 2(F), 3(K)

MVF_N_vars = 70;  % Number of total variable coulmns to allow 


%% User Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
site = 'SNQ15';
Forcing_N_in  = 'OctMar';
Forcing_N_out = 'OctMar';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('QCing for %s \n',Forcing_N_out);
file_in = strcat(site,'_RAW.mat');
file_out = strcat(site,'_QC.mat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Directories
homedir = 'G:\ncar\d1\data\Point_Snow_Sites\'; % REPLACE YOUR HOME DIR HERE
rawdir = fullfile(homedir,'Raw',site,Forcing_N_in);
QCdir = fullfile(homedir,'QC',site,Forcing_N_out);
if exist(QCdir,'dir') ~= 7
    mkdir(QCdir)
end
cd(QCdir)
%% Load Raw MVF data
load(strcat(rawdir,'\',file_in)); 

%% Quality Control Input Data
Data_QC = NaN(size(Data_Raw)); % QC'ed data (only includes variables QC'ed below)

% Find columns with all NaN
noData = find(sum(isnan(Data_Raw),1) == size(Data_Raw,1));

%% MVF - 1 - Air Temperature A

disp('Air Temperature A')
VAR_N = 1; 
Var_file_save = '\QC_Temp.mat';

% Check we have any data
if ~any(noData == VAR_N) % There is some data for this variable
    if exist(strcat(QCdir,Var_file_save)) == 2 %#ok<*EXIST>
        disp('Previous Quality Controled Temp A found')
        u_inp = user_inp(cellstr('Enter 1 - Use this data. 2 - Create new QC data.'), 1, 2);
    else
        u_inp = 2;
    end

    if u_inp == 2
        % Get Raw data
        Var_raw = Data_Raw(:,VAR_N);
        % QC (Interactive)
        [QC_Temp, accepted] = data_QC_NW_edits(Var_raw, 1, TIME_out, 1);
        % Save Var QC data
        eval(['save ' strcat(QCdir,Var_file_save) ' QC_Temp']);
    elseif u_inp == 1
       load(strcat(QCdir,Var_file_save)) 
    end
    Data_QC(:,VAR_N) = QC_Temp;
else % There is no data for this variable
    Data_QC(:,VAR_N) = NaN; 
end
clc;

%% MVF - 31 - Air Temperature B

disp('Air Temperature B')
VAR_N = 31; 
Var_file_save = '\QC_TempB.mat';

% Check we have any data
if ~any(noData == VAR_N) % There is some data for this variable
    if exist(strcat(QCdir,Var_file_save)) == 2
        disp('Previous Quality Controled Temp B found')
        u_inp = user_inp(cellstr('Enter 1 - Use this data. 2 - Create new QC data.'), 1, 2);
    else
        u_inp = 2;
    end

    if u_inp == 2
        % Get Raw data
        Var_raw = Data_Raw(:,VAR_N);
        % QC (Interactive)
        [QC_TempB, accepted] = data_QC_NW_edits(Var_raw, 1, TIME_out, 1);
        % Save Var QC data
        eval(['save ' strcat(QCdir,Var_file_save) ' QC_TempB']);
    elseif u_inp == 1
       load(strcat(QCdir,Var_file_save)) 
    end
    Data_QC(:,VAR_N) = QC_TempB;
else % There is no data for this variable
    Data_QC(:,VAR_N) = NaN; 
end
clc;

%% MVF - 42 - Air Temperature C

disp('Air Temperature C')
VAR_N = 42; 
Var_file_save = '\QC_TempC.mat';

% Check we have any data
if ~any(noData == VAR_N) % There is some data for this variable
    if exist(strcat(QCdir,Var_file_save)) == 2
        disp('Previous Quality Controled Temp C found')
        u_inp = user_inp(cellstr('Enter 1 - Use this data. 2 - Create new QC data.'), 1, 2);
    else
        u_inp = 2;
    end

    if u_inp == 2
        % Get Raw data
        Var_raw = Data_Raw(:,VAR_N);
        % QC (Interactive)
        [QC_TempC, accepted] = data_QC_NW_edits(Var_raw, 1, TIME_out, 1);
        % Save Var QC data
        eval(['save ' strcat(QCdir,Var_file_save) ' QC_TempC']);
    elseif u_inp == 1
       load(strcat(QCdir,Var_file_save)) 
    end
    Data_QC(:,VAR_N) = QC_TempC;
else % There is no data for this variable
    Data_QC(:,VAR_N) = NaN; 
end
clc;

%% MVF - 43 - Air Temperature D

disp('Air Temperature D')
VAR_N = 43; 
Var_file_save = '\QC_TempD.mat';

% Check we have any data
if ~any(noData == VAR_N) % There is some data for this variable
    if exist(strcat(QCdir,Var_file_save)) == 2
        disp('Previous Quality Controled Temp D found')
        u_inp = user_inp(cellstr('Enter 1 - Use this data. 2 - Create new QC data.'), 1, 2);
    else
        u_inp = 2;
    end

    if u_inp == 2
        % Get Raw data
        Var_raw = Data_Raw(:,VAR_N);
        % QC (Interactive)
        [QC_TempD, accepted] = data_QC_NW_edits(Var_raw, 1, TIME_out, 1);
        % Save Var QC data
        eval(['save ' strcat(QCdir,Var_file_save) ' QC_TempD']);
    elseif u_inp == 1
       load(strcat(QCdir,Var_file_save)) 
    end
    Data_QC(:,VAR_N) = QC_TempD;
else % There is no data for this variable
    Data_QC(:,VAR_N) = NaN; 
end
clc;


%% MVF - 17 - Dew Point

disp('Dew Point')
VAR_N = 17; 
Var_file_save = '\QC_DP.mat';

% Check we have any data
if ~any(noData == VAR_N) % There is some data for this variable
    if exist(strcat(QCdir,Var_file_save)) == 2
        disp('Previous Quality Controled Dew Point found')
        u_inp = user_inp(cellstr('Enter 1 - Use this data. 2 - Create new QC data.'), 1, 2);
    else
        u_inp = 2;
    end

    if u_inp == 2
        % Get Raw data
        Var_raw = Data_Raw(:,VAR_N);
        % QC (Interactive)
        [QC_DP, accepted] = data_QC_NW_edits(Var_raw, 1, TIME_out, 1);
        % Save Var QC data
        eval(['save ' strcat(QCdir,Var_file_save) ' QC_DP']);
    elseif u_inp == 1
       load(strcat(QCdir,Var_file_save)) 
    end
    Data_QC(:,VAR_N) = QC_DP;
else % There is no data for this variable
    Data_QC(:,VAR_N) = NaN; 
end
clc;

%% MVF - 2 - Precipitation
disp('Precipitation')
VAR_N = 2;
Var_file_save = '\QC_Precip.mat';

% Check we have any data
if ~any(noData == VAR_N) % There is some data for this variable
    if exist(strcat(QCdir,Var_file_save)) == 2
        disp('Previous Quality Controled Precip found')
        u_inp = user_inp(cellstr('Enter 1 - Use this data. 2 - Create new QC data.'), 1, 2);
    else
        u_inp = 2;
    end

    if u_inp == 2
        % Get Raw data
        Var_raw = Data_Raw(:,VAR_N);
        % data_QC wants mm not m
        Var_raw = Var_raw .* 1000; 
        % QC (Interactive)
        [QC_Precip, accepted] = data_QC_NW_edits(Var_raw, 1, TIME_out, 2);
        QC_Precip = QC_Precip ./ 1000; % mm to m (MVF UNITS)
        % Save Var QC data
        eval(['save ' strcat(QCdir,Var_file_save) ' QC_Precip']);
    elseif u_inp == 1
       load(strcat(QCdir,Var_file_save)) 
    end
    Data_QC(:,VAR_N) = QC_Precip;
else % There is no data for this variable
    Data_QC(:,VAR_N) = NaN; 
end 
clc;

%% MVF - 55 - Incremental Precipitation B
disp('Incremental Precipitation B')
VAR_N = 55;
Var_file_save = '\QC_Precip_B.mat';

% Check we have any data
if ~any(noData == VAR_N) % There is some data for this variable
    if exist(strcat(QCdir,Var_file_save)) == 2
        disp('Previous Quality Controled Incremental Precipitation B found')
        u_inp = user_inp(cellstr('Enter 1 - Use this data. 2 - Create new QC data.'), 1, 2);
    else
        u_inp = 2;
    end

    if u_inp == 2
        % Get Raw data
        Var_raw = Data_Raw(:,VAR_N);
        % data_QC wants mm not m
        Var_raw = Var_raw .* 1000; 
        % QC (Interactive)
        [QC_Precip_B, accepted] = data_QC_NW_edits(Var_raw, 1, TIME_out, 2);
        QC_Precip_B = QC_Precip_B ./ 1000; % mm to m (MVF UNITS)
        % Save Var QC data
        eval(['save ' strcat(QCdir,Var_file_save) ' QC_Precip_B']);
    elseif u_inp == 1
       load(strcat(QCdir,Var_file_save)) 
    end
    Data_QC(:,VAR_N) = QC_Precip_B;
else % There is no data for this variable
    Data_QC(:,VAR_N) = NaN; 
end 
clc;

%% MVF - 3 - Precipitation
disp('Precipitation')
VAR_N = 3;
VAR_N_out = 2; % to incremental
Var_file_save = '\QC_Precip.mat';

% Check we have any data
if ~any(noData == VAR_N) % There is some data for this variable
    if exist(strcat(QCdir,Var_file_save)) == 2
        disp('Previous Quality Controled accumluated found')
        u_inp = user_inp(cellstr('Enter 1 - Use this data. 2 - Create new QC data.'), 1, 2);
    else
        u_inp = 2;
    end

    if u_inp == 2
        % Get Raw data
        Var_raw = Data_Raw(:,VAR_N);
        % data_QC wants mm not m
        Var_raw = Var_raw .* 1000; 
        
        disp('WARNING: simple conversion form accumlative to incremental precipitation, check it!')
        disp('Max precip amount hardcoded at ')      
        
        % Convert from accumulative to incremental
        Var_raw_diff = [0; nandiff(Var_raw)];
        Var_raw_diff(Var_raw_diff<0) = NaN; % Remove negative accumulations
        Var_raw_diff(Var_raw_diff>100) = NaN; % Remove unrealistic accumulations
        
        % QC (Interactive)
        [QC_Precip, accepted] = data_QC_NW_edits(Var_raw, 1, TIME_out, 2);
        QC_Precip = QC_Precip ./ 1000; % mm to m (MVF UNITS)
        % Save Var QC data
        eval(['save ' strcat(QCdir,Var_file_save) ' QC_Precip']);
    elseif u_inp == 1
       load(strcat(QCdir,Var_file_save)) 
    end
    Data_QC(:,VAR_N_out) = QC_Precip; % Saving out to incremental
else % There is no data for this variable
    Data_QC(:,VAR_N) = NaN; 
end 
clc;

%% MVF - 4 - Relative Humidity
disp('Relative Humidity A')
VAR_N = 4;
Var_file_save = '\QC_RH.mat';

% Check we have any data
if ~any(noData == VAR_N) % There is some data for this variable
    if exist(strcat(QCdir,Var_file_save)) == 2
        disp('Previous Quality Controled Relative Humidity A')
        u_inp = user_inp(cellstr('Enter 1 - Use this data. 2 - Create new QC data.'), 1, 2);
    else
        u_inp = 2;
    end

    if u_inp == 2
        % Get Raw data
        Var_raw = Data_Raw(:,VAR_N);
        % QC (Interactive)
        [QC_RH, accepted] = data_QC_NW_edits(Var_raw, 1, TIME_out, 3);
        % Save Var QC data
        eval(['save ' strcat(QCdir,Var_file_save) ' QC_RH']);
    elseif u_inp == 1
       load(strcat(QCdir,Var_file_save)) 
    end
    Data_QC(:,VAR_N) = QC_RH;
else % There is no data for this variable
    Data_QC(:,VAR_N) = NaN; 
end 
clc;

%% MVF - 32 - Relative Humidity
disp('Relative Humidity B')
VAR_N = 32;
Var_file_save = '\QC_RHB.mat';

% Check we have any data
if ~any(noData == VAR_N) % There is some data for this variable
    if exist(strcat(QCdir,Var_file_save)) == 2
        disp('Previous Quality Controled Relative Humidity B')
        u_inp = user_inp(cellstr('Enter 1 - Use this data. 2 - Create new QC data.'), 1, 2);
    else
        u_inp = 2;
    end

    if u_inp == 2
        % Get Raw data
        Var_raw = Data_Raw(:,VAR_N);
        % QC (Interactive)
        [QC_RHB, accepted] = data_QC_NW_edits(Var_raw, 1, TIME_out, 3);
        % Save Var QC data
        eval(['save ' strcat(QCdir,Var_file_save) ' QC_RHB']);
    elseif u_inp == 1
       load(strcat(QCdir,Var_file_save)) 
    end
    Data_QC(:,VAR_N) = QC_RHB;
else % There is no data for this variable
    Data_QC(:,VAR_N) = NaN; 
end 
clc;

%% MVF - 5 - Wind Speed A
disp('Wind Speed A')
VAR_N = 5;
Var_file_save = '\QC_WSA.mat';

% Check we have any data
if ~any(noData == VAR_N) % There is some data for this variable
    if exist(strcat(QCdir,Var_file_save)) == 2
        disp('Previous Quality Controled Wind Speed A found')
        u_inp = user_inp(cellstr('Enter 1 - Use this data. 2 - Create new QC data.'), 1, 2);
    else
        u_inp = 2;
    end

    if u_inp == 2
        % Get Raw data
        Var_raw = Data_Raw(:,VAR_N);
        % QC (Interactive)
        [QC_WSA, accepted] = data_QC_NW_edits(Var_raw, 1, TIME_out, 4);
        if accepted == 0
            disp('Exiting because QC was aborted')
            return
        end
        % Save Var QC data
        eval(['save ' strcat(QCdir,Var_file_save) ' QC_WSA']);
    elseif u_inp == 1
       load(strcat(QCdir,Var_file_save)) 
    end
    Data_QC(:,VAR_N) = QC_WSA;
else % There is no data for this variable
    Data_QC(:,VAR_N) = NaN; 
end 
clc;

%% MVF - 6 - Wind Speed B
disp('Wind Speed B')
VAR_N = 6;
Var_file_save = '\QC_WSB.mat';

% Check we have any data
if ~any(noData == VAR_N) % There is some data for this variable
    if exist(strcat(QCdir,Var_file_save)) == 2
        disp('Previous Quality Controled Wind Speed B found')
        u_inp = user_inp(cellstr('Enter 1 - Use this data. 2 - Create new QC data.'), 1, 2);
    else
        u_inp = 2;
    end

    if u_inp == 2
        % Get Raw data
        Var_raw = Data_Raw(:,VAR_N);
        % QC (Interactive)
        [QC_WSB, accepted] = data_QC_NW_edits(Var_raw, 1, TIME_out, 4);
        if accepted == 0
            disp('Exiting because QC was aborted')
            return
        end
        % Save Var QC data
        eval(['save ' strcat(QCdir,Var_file_save) ' QC_WSB']);
    elseif u_inp == 1
       load(strcat(QCdir,Var_file_save)) 
    end
    Data_QC(:,VAR_N) = QC_WSB;
else % There is no data for this variable
    Data_QC(:,VAR_N) = NaN; 
end 
clc;

%% MVF - 14 - Pressure adjusted to sea-level
disp('Adjusted Pressure')
VAR_N = 14;
VAR_N_alt = 16; % Where we put readjusted pressure
Var_file_save = '\QC_Press.mat';

% Check we have any data
if ~any(noData == VAR_N) % There is some data for this variable
    if exist(strcat(QCdir,Var_file_save)) == 2
        disp('Previous Quality Controled Pressure found')
        u_inp = user_inp(cellstr('Enter 1 - Use this data. 2 - Create new QC data.'), 1, 2);
    else
        u_inp = 2;
    end

    if u_inp == 2
        % Get Raw data
        Var_raw = Data_Raw(:,VAR_N); 
        % QC (Interactive)
        [QC_Press, accepted] = data_QC_NW_edits(Var_raw, 1, TIME_out, 7);
        if accepted == 0
            disp('Exiting because QC was aborted')
            return
        end
        % Save Var QC data
        eval(['save ' strcat(QCdir,Var_file_save) ' QC_Press']);
    elseif u_inp == 1
       load(strcat(QCdir,Var_file_save)) 
    end
    Data_QC(:,VAR_N) = QC_Press;
else % There is no data for this variable
    Data_QC(:,VAR_N) = NaN; 
end     

if ~any(noData == VAR_N)
    disp('Would you like to convert Adjusted Pressure to Unadjusted Pressure B? (will not overwrite any data)')
    u_inpp = user_inp(cellstr('Enter 1 - Yes. 2 - No.'), 1, 2);
    if u_inpp == 1
        % Check where to put new converted Pressure in MVF
        if sum(~isnan(Data_Raw(:,VAR_N_alt))) > 0 % then there is already data here
            error('There is already raw data in VAR_N_alt column of MVF, please choose another')
        else
           % Convert back from adjusted to unadjusted pressure
           Data_QC(:,VAR_N_alt) = Press_from_sea_2_level(Data_QC(:,VAR_N),site_E);
           fprintf('Converting sea-level presure to %3.2fm-level pressure (unadjusted) \n',site_E); 
        end
    end
end
clc;

%% MVF - 46 - Pressure unadjusted
disp('Unadjusted Pressure')
VAR_N = 46;
Var_file_save = '\QC_PressUnAdj.mat';
% Check we have any data
if ~any(noData == VAR_N) % There is some data for this variable
    if exist(strcat(QCdir,Var_file_save)) == 2
        disp('Previous Quality Controled Unadjusted Pressure found')
        u_inp = user_inp(cellstr('Enter 1 - Use this data. 2 - Create new QC data. 3 - Continue with previous QC data'), 1, 2, 3);
    else
        u_inp = 2;
    end

    if u_inp == 2 || u_inp == 3
        if u_inp == 2
            % Get Raw data
            Var_raw = Data_Raw(:,VAR_N); 
        elseif u_inp == 3 % Get previous Filled data
            load(strcat(QCdir,Var_file_save)) 
            Var_raw = QC_PressUnAdj;
        end
            
        % QC (Interactive)
        [QC_PressUnAdj, accepted] = data_QC_NW_edits(Var_raw, 1, TIME_out, 7);
        if accepted == 0
            disp('Exiting because QC was aborted')
            return
        end
        % Save Var QC data
        eval(['save ' strcat(QCdir,Var_file_save) ' QC_PressUnAdj']);
    elseif u_inp == 1
       load(strcat(QCdir,Var_file_save)) 
    end
    Data_QC(:,VAR_N) = QC_PressUnAdj;
else % There is no data for this variable
    Data_QC(:,VAR_N) = NaN; 
end 
clc;

%% MVF - 9 - Shortwave downward
disp('Shortwave downward')
VAR_N = 9;
Var_file_save = '\QC_SWd.mat';

% Check we have any data
if ~any(noData == VAR_N) % There is some data for this variable
    if exist(strcat(QCdir,Var_file_save)) == 2
        disp('Previous Quality Controled SW down found')
        u_inp = user_inp(cellstr('Enter 1 - Use this data. 2 - Create new QC data.'), 1, 2);
    else
        u_inp = 2;
    end

    if u_inp == 2
        % Get Raw data
        Var_raw = Data_Raw(:,VAR_N); 
        % QC (Interactive)
        [QC_SWd, accepted] = data_QC_NW_edits(Var_raw, 1, TIME_out, 5);
        if accepted == 0
            disp('Exiting because QC was aborted')
            return
        end
        % Save Var QC data
        eval(['save ' strcat(QCdir,Var_file_save) ' QC_SWd']);
    elseif u_inp == 1
       load(strcat(QCdir,Var_file_save)) 
    end
    Data_QC(:,VAR_N) = QC_SWd;
else % There is no data for this variable
    Data_QC(:,VAR_N) = NaN; 
end
clc;


%% MVF - 11 - Longwave downward
disp('Longwave downward')
VAR_N = 11;
Var_file_save = '\QC_LWd.mat';

% Check we have any data
if ~any(noData == VAR_N) % There is some data for this variable
    if exist(strcat(QCdir,Var_file_save)) == 2
        disp('Previous Quality Controled LW down found')
        u_inp = user_inp(cellstr('Enter 1 - Use this data. 2 - Create new QC data.'), 1, 2);
    else
        u_inp = 2;
    end

    if u_inp == 2
        % Get Raw data
        Var_raw = Data_Raw(:,VAR_N); 
        % QC (Interactive)
        [QC_LWd, accepted] = data_QC_NW_edits(Var_raw, 1, TIME_out, 8);
        if accepted == 0
            disp('Exiting because QC was aborted')
            return
        end
        % Save Var QC data
        eval(['save ' strcat(QCdir,Var_file_save) ' QC_LWd']);
    elseif u_inp == 1
       load(strcat(QCdir,Var_file_save)) 
    end
    Data_QC(:,VAR_N) = QC_LWd;
else % There is no data for this variable
    Data_QC(:,VAR_N) = NaN; 
end
clc;


%% END FORCING




%% START VERIFICATION

%% MVF - 10 - Shortwave upward
disp('Shortwave upward')
VAR_N = 10;
Var_file_save = '\QC_SWu.mat';

% Check we have any data
if ~any(noData == VAR_N) % There is some data for this variable
    if exist(strcat(QCdir,Var_file_save)) == 2
        disp('Previous Quality Controled SW up found')
        u_inp = user_inp(cellstr('Enter 1 - Use this data. 2 - Create new QC data.'), 1, 2);
    else
        u_inp = 2;
    end

    if u_inp == 2
        % Get Raw data
        Var_raw = Data_Raw(:,VAR_N); 
        % QC (Interactive)
        [QC_SWu, accepted] = data_QC_NW_edits(Var_raw, 1, TIME_out, 5);
        if accepted == 0
            disp('Exiting because QC was aborted')
            return
        end
        % Save Var QC data
        eval(['save ' strcat(QCdir,Var_file_save) ' QC_SWu']);
    elseif u_inp == 1
       load(strcat(QCdir,Var_file_save)) 
    end
    Data_QC(:,VAR_N) = QC_SWu;
else % There is no data for this variable
    Data_QC(:,VAR_N) = NaN; 
end
clc;


%% MVF - 12 - Longwave upward
disp('Longwave upward')
VAR_N = 12;
Var_file_save = '\QC_LWu.mat';

% Check we have any data
if ~any(noData == VAR_N) % There is some data for this variable
    if exist(strcat(QCdir,Var_file_save)) == 2
        disp('Previous Quality Controled LW up found')
        u_inp = user_inp(cellstr('Enter 1 - Use this data. 2 - Create new QC data.'), 1, 2);
    else
        u_inp = 2;
    end

    if u_inp == 2
        % Get Raw data
        Var_raw = Data_Raw(:,VAR_N); 
        % QC (Interactive)
        [QC_LWu, accepted] = data_QC_NW_edits(Var_raw, 1, TIME_out, 8);
        if accepted == 0
            disp('Exiting because QC was aborted')
            return
        end
        % Save Var QC data
        eval(['save ' strcat(QCdir,Var_file_save) ' QC_LWu']);
    elseif u_inp == 1
       load(strcat(QCdir,Var_file_save)) 
    end
    Data_QC(:,VAR_N) = QC_LWu;
else % There is no data for this variable
    Data_QC(:,VAR_N) = NaN; 
end    
clc;

%% MVF - 25 - Sensible Heat Flux 
disp('Sensible Heat Flux ')
VAR_N = 25;
Var_file_save = '\QC_SH.mat';

% Check we have any data
if ~any(noData == VAR_N) % There is some data for this variable
    if exist(strcat(QCdir,Var_file_save)) == 2
        disp('Previous Quality Controled Sensible Heat Flux down found')
        u_inp = user_inp(cellstr('Enter 1 - Use this data. 2 - Create new QC data.'), 1, 2);
    else
        u_inp = 2;
    end

    if u_inp == 2
        % Get Raw data
        Var_raw = Data_Raw(:,VAR_N); 
        % QC (Interactive)
        [QC_SH, accepted] = data_QC_NW_edits(Var_raw, 1, TIME_out, 5);
        if accepted == 0
            disp('Exiting because QC was aborted')
            return
        end
        % Save Var QC data
        eval(['save ' strcat(QCdir,Var_file_save) ' QC_SH']);
    elseif u_inp == 1
       load(strcat(QCdir,Var_file_save)) 
    end
    Data_QC(:,VAR_N) = QC_SH;
else % There is no data for this variable
    Data_QC(:,VAR_N) = NaN; 
end      
clc;

%% MVF - 26 - Latent Heat Flux 
disp('Latent Heat Flux ')
VAR_N = 26;
Var_file_save = '\QC_LH.mat';

% Check we have any data
if ~any(noData == VAR_N) % There is some data for this variable
    if exist(strcat(QCdir,Var_file_save)) == 2
        disp('Previous Quality Controled Latent Heat Flux down found')
        u_inp = user_inp(cellstr('Enter 1 - Use this data. 2 - Create new QC data.'), 1, 2);
    else
        u_inp = 2;
    end

    if u_inp == 2
        % Get Raw data
        Var_raw = Data_Raw(:,VAR_N); 
        % QC (Interactive)
        [QC_LH, accepted] = data_QC_NW_edits(Var_raw, 1, TIME_out, 5);
        if accepted == 0
            disp('Exiting because QC was aborted')
            return
        end
        % Save Var QC data
        eval(['save ' strcat(QCdir,Var_file_save) ' QC_LH']);
    elseif u_inp == 1
       load(strcat(QCdir,Var_file_save)) 
    end
    Data_QC(:,VAR_N) = QC_LH;
else % There is no data for this variable
    Data_QC(:,VAR_N) = NaN; 
end      
clc;

%% MVF - 33 - Soil Heat Flux
disp('Soil Heat Flux ')
VAR_N = 33;
Var_file_save = '\SoilHFavg.mat';

% Check we have any data
if ~any(noData == VAR_N) % There is some data for this variable
    if exist(strcat(QCdir,Var_file_save)) == 2
        disp('Previous Quality Controled Soil Heat Flux found')
        u_inp = user_inp(cellstr('Enter 1 - Use this data. 2 - Create new QC data.'), 1, 2);
    else
        u_inp = 2;
    end

    if u_inp == 2
        % Get Raw data
        Var_raw = Data_Raw(:,VAR_N); 
        % QC (Interactive)
        [SoilHFavg, accepted] = data_QC_NW_edits(Var_raw, 1, TIME_out, 5);
        if accepted == 0
            disp('Exiting because QC was aborted')
            return
        end
        % Save Var QC data
        eval(['save ' strcat(QCdir,Var_file_save) ' SoilHFavg']);
    elseif u_inp == 1
       load(strcat(QCdir,Var_file_save)) 
    end
    Data_QC(:,VAR_N) = SoilHFavg;
else % There is no data for this variable
    Data_QC(:,VAR_N) = NaN; 
end      
clc;


%% MVF - 22 - Surface Temperature
disp('Surface Temperature A')
VAR_N = 22;
Var_file_save = '\QC_ST.mat';

% Check we have any data
if ~any(noData == VAR_N) % There is some data for this variable
    if exist(strcat(QCdir,Var_file_save)) == 2
        disp('Previous Quality Controled Surface Temperature A down found')
        u_inp = user_inp(cellstr('Enter 1 - Use this data. 2 - Create new QC data.'), 1, 2);
    else
        u_inp = 2;
    end

    if u_inp == 2
        % Get Raw data
        Var_raw = Data_Raw(:,VAR_N); 
        % QC (Interactive)
        [QC_ST, accepted] = data_QC_NW_edits(Var_raw, 1, TIME_out, 1);
        if accepted == 0
            disp('Exiting because QC was aborted')
            return
        end
        % Save Var QC data
        eval(['save ' strcat(QCdir,Var_file_save) ' QC_ST']);
    elseif u_inp == 1
       load(strcat(QCdir,Var_file_save)) 
    end
    Data_QC(:,VAR_N) = QC_ST;
else % There is no data for this variable
    Data_QC(:,VAR_N) = NaN; 
end      
clc;

%% MVF - 51 - Surface Temperature B 
disp('Surface Temperature B ')
VAR_N = 51;
Var_file_save = '\QC_ST_B.mat';

% Check we have any data
if ~any(noData == VAR_N) % There is some data for this variable
    if exist(strcat(QCdir,Var_file_save)) == 2
        disp('Previous Quality Controled Surface Temperature B  down found')
        u_inp = user_inp(cellstr('Enter 1 - Use this data. 2 - Create new QC data.'), 1, 2);
    else
        u_inp = 2;
    end

    if u_inp == 2
        % Get Raw data
        Var_raw = Data_Raw(:,VAR_N); 
        % QC (Interactive)
        [QC_ST_B, accepted] = data_QC_NW_edits(Var_raw, 1, TIME_out, 1);
        if accepted == 0
            disp('Exiting because QC was aborted')
            return
        end
        % Save Var QC data
        eval(['save ' strcat(QCdir,Var_file_save) ' QC_ST_B']);
    elseif u_inp == 1
       load(strcat(QCdir,Var_file_save)) 
    end
    Data_QC(:,VAR_N) = QC_ST_B;
else % There is no data for this variable
    Data_QC(:,VAR_N) = NaN; 
end      
clc;

%% MVF - 62 - Surface Temperature C 
disp('Surface Temperature C ')
VAR_N = 62;
Var_file_save = '\QC_ST_C.mat';

% Check we have any data
if ~any(noData == VAR_N) % There is some data for this variable
    if exist(strcat(QCdir,Var_file_save)) == 2
        disp('Previous Quality Controled Surface Temperature C  down found')
        u_inp = user_inp(cellstr('Enter 1 - Use this data. 2 - Create new QC data.'), 1, 2);
    else
        u_inp = 2;
    end

    if u_inp == 2
        % Get Raw data
        Var_raw = Data_Raw(:,VAR_N); 
        % QC (Interactive)
        [QC_ST_C, accepted] = data_QC_NW_edits(Var_raw, 1, TIME_out, 1);
        if accepted == 0
            disp('Exiting because QC was aborted')
            return
        end
        % Save Var QC data
        eval(['save ' strcat(QCdir,Var_file_save) ' QC_ST_C']);
    elseif u_inp == 1
       load(strcat(QCdir,Var_file_save)) 
    end
    Data_QC(:,VAR_N) = QC_ST_C;
else % There is no data for this variable
    Data_QC(:,VAR_N) = NaN; 
end      
clc;

%% MVF - 37 - Soil Temperature A
disp('Soil Temperature A ')
VAR_N = 37;
Var_file_save = '\SoilTA.mat';

% Check we have any data
if ~any(noData == VAR_N) % There is some data for this variable
    if exist(strcat(QCdir,Var_file_save)) == 2
        disp('Previous Quality Controled Soil Temperature A  down found')
        u_inp = user_inp(cellstr('Enter 1 - Use this data. 2 - Create new QC data.'), 1, 2);
    else
        u_inp = 2;
    end

    if u_inp == 2
        % Get Raw data
        Var_raw = Data_Raw(:,VAR_N); 
        % QC (Interactive)
        [SoilTA, accepted] = data_QC_NW_edits(Var_raw, 1, TIME_out, 1);
        if accepted == 0
            disp('Exiting because QC was aborted')
            return
        end
        % Save Var QC data
        eval(['save ' strcat(QCdir,Var_file_save) ' SoilTA']);
    elseif u_inp == 1
       load(strcat(QCdir,Var_file_save)) 
    end
    Data_QC(:,VAR_N) = SoilTA;
else % There is no data for this variable
    Data_QC(:,VAR_N) = NaN; 
end      
clc;

%% MVF - 38 - Soil Temperature B
disp('Soil Temperature B ')
VAR_N = 38;
Var_file_save = '\SoilTB.mat';

% Check we have any data
if ~any(noData == VAR_N) % There is some data for this variable
    if exist(strcat(QCdir,Var_file_save)) == 2
        disp('Previous Quality Controled Soil Temperature B  down found')
        u_inp = user_inp(cellstr('Enter 1 - Use this data. 2 - Create new QC data.'), 1, 2);
    else
        u_inp = 2;
    end

    if u_inp == 2
        % Get Raw data
        Var_raw = Data_Raw(:,VAR_N); 
        % QC (Interactive)
        [SoilTB, accepted] = data_QC_NW_edits(Var_raw, 1, TIME_out, 1);
        if accepted == 0
            disp('Exiting because QC was aborted')
            return
        end
        % Save Var QC data
        eval(['save ' strcat(QCdir,Var_file_save) ' SoilTB']);
    elseif u_inp == 1
       load(strcat(QCdir,Var_file_save)) 
    end
    Data_QC(:,VAR_N) = SoilTB;
else % There is no data for this variable
    Data_QC(:,VAR_N) = NaN; 
end      
clc;

%% MVF - 39 - Soil Temperature C
disp('Soil Temperature C ')
VAR_N = 39;
Var_file_save = '\SoilTC.mat';

% Check we have any data
if ~any(noData == VAR_N) % There is some data for this variable
    if exist(strcat(QCdir,Var_file_save)) == 2
        disp('Previous Quality Controled Soil Temperature C  down found')
        u_inp = user_inp(cellstr('Enter 1 - Use this data. 2 - Create new QC data.'), 1, 2);
    else
        u_inp = 2;
    end

    if u_inp == 2
        % Get Raw data
        Var_raw = Data_Raw(:,VAR_N); 
        % QC (Interactive)
        [SoilTC, accepted] = data_QC_NW_edits(Var_raw, 1, TIME_out, 1);
        if accepted == 0
            disp('Exiting because QC was aborted')
            return
        end
        % Save Var QC data
        eval(['save ' strcat(QCdir,Var_file_save) ' SoilTC']);
    elseif u_inp == 1
       load(strcat(QCdir,Var_file_save)) 
    end
    Data_QC(:,VAR_N) = SoilTC;
else % There is no data for this variable
    Data_QC(:,VAR_N) = NaN; 
end      
clc;

%% MVF - 40 - Soil Temperature D
disp('Soil Temperature D ')
VAR_N = 40;
Var_file_save = '\SoilTD.mat';

% Check we have any data
if ~any(noData == VAR_N) % There is some data for this variable
    if exist(strcat(QCdir,Var_file_save)) == 2
        disp('Previous Quality Controled Soil Temperature D  down found')
        u_inp = user_inp(cellstr('Enter 1 - Use this data. 2 - Create new QC data.'), 1, 2);
    else
        u_inp = 2;
    end

    if u_inp == 2
        % Get Raw data
        Var_raw = Data_Raw(:,VAR_N); 
        % QC (Interactive)
        [SoilTD, accepted] = data_QC_NW_edits(Var_raw, 1, TIME_out, 1);
        if accepted == 0
            disp('Exiting because QC was aborted')
            return
        end
        % Save Var QC data
        eval(['save ' strcat(QCdir,Var_file_save) ' SoilTD']);
    elseif u_inp == 1
       load(strcat(QCdir,Var_file_save)) 
    end
    Data_QC(:,VAR_N) = SoilTD;
else % There is no data for this variable
    Data_QC(:,VAR_N) = NaN; 
end      
clc;

%% MVF - 41 - Soil Temperature E
disp('Soil Temperature E ')
VAR_N = 41;
Var_file_save = '\SoilTE.mat';

% Check we have any data
if ~any(noData == VAR_N) % There is some data for this variable
    if exist(strcat(QCdir,Var_file_save)) == 2
        disp('Previous Quality Controled Soil Temperature E  down found')
        u_inp = user_inp(cellstr('Enter 1 - Use this data. 2 - Create new QC data.'), 1, 2);
    else
        u_inp = 2;
    end

    if u_inp == 2
        % Get Raw data
        Var_raw = Data_Raw(:,VAR_N); 
        % QC (Interactive)
        [SoilTE, accepted] = data_QC_NW_edits(Var_raw, 1, TIME_out, 1);
        if accepted == 0
            disp('Exiting because QC was aborted')
            return
        end
        % Save Var QC data
        eval(['save ' strcat(QCdir,Var_file_save) ' SoilTE']);
    elseif u_inp == 1
       load(strcat(QCdir,Var_file_save)) 
    end
    Data_QC(:,VAR_N) = SoilTE;
else % There is no data for this variable
    Data_QC(:,VAR_N) = NaN; 
end      
clc;


%% MVF - 56 - Soil Moisture A 
disp('Soil Moisture A ')
VAR_N = 56;
Var_file_save = '\SoilMA.mat';

% Check we have any data
if ~any(noData == VAR_N) % There is some data for this variable
    if exist(strcat(QCdir,Var_file_save)) == 2
        disp('Previous Quality Controled Soil Moisture A   down found')
        u_inp = user_inp(cellstr('Enter 1 - Use this data. 2 - Create new QC data.'), 1, 2);
    else
        u_inp = 2;
    end

    if u_inp == 2
        % Get Raw data
        Var_raw = Data_Raw(:,VAR_N); 
        % QC (Interactive)
        [SoilMA, accepted] = data_QC_NW_edits(Var_raw, 1, TIME_out, 10);
        if accepted == 0
            disp('Exiting because QC was aborted')
            return
        end
        % Save Var QC data
        eval(['save ' strcat(QCdir,Var_file_save) ' SoilMA']);
    elseif u_inp == 1
       load(strcat(QCdir,Var_file_save)) 
    end
    Data_QC(:,VAR_N) = SoilMA;
else % There is no data for this variable
    Data_QC(:,VAR_N) = NaN; 
end      
clc;

%% MVF - 57 - Soil Moisture B 
disp('Soil Moisture B ')
VAR_N = 57;
Var_file_save = '\SoilMB.mat';

% Check we have any data
if ~any(noData == VAR_N) % There is some data for this variable
    if exist(strcat(QCdir,Var_file_save)) == 2
        disp('Previous Quality Controled Soil Moisture B   down found')
        u_inp = user_inp(cellstr('Enter 1 - Use this data. 2 - Create new QC data.'), 1, 2);
    else
        u_inp = 2;
    end

    if u_inp == 2
        % Get Raw data
        Var_raw = Data_Raw(:,VAR_N); 
        % QC (Interactive)
        [SoilMB, accepted] = data_QC_NW_edits(Var_raw, 1, TIME_out, 10);
        if accepted == 0
            disp('Exiting because QC was aborted')
            return
        end
        % Save Var QC data
        eval(['save ' strcat(QCdir,Var_file_save) ' SoilMB']);
    elseif u_inp == 1
       load(strcat(QCdir,Var_file_save)) 
    end
    Data_QC(:,VAR_N) = SoilMB;
else % There is no data for this variable
    Data_QC(:,VAR_N) = NaN; 
end      
clc;

%% MVF - 58 - Soil Moisture C 
disp('Soil Moisture C ')
VAR_N = 58;
Var_file_save = '\SoilMC.mat';

% Check we have any data
if ~any(noData == VAR_N) % There is some data for this variable
    if exist(strcat(QCdir,Var_file_save)) == 2
        disp('Previous Quality Controled Soil Moisture C   down found')
        u_inp = user_inp(cellstr('Enter 1 - Use this data. 2 - Create new QC data.'), 1, 2);
    else
        u_inp = 2;
    end

    if u_inp == 2
        % Get Raw data
        Var_raw = Data_Raw(:,VAR_N); 
        % QC (Interactive)
        [SoilMC, accepted] = data_QC_NW_edits(Var_raw, 1, TIME_out, 10);
        if accepted == 0
            disp('Exiting because QC was aborted')
            return
        end
        % Save Var QC data
        eval(['save ' strcat(QCdir,Var_file_save) ' SoilMC']);
    elseif u_inp == 1
       load(strcat(QCdir,Var_file_save)) 
    end
    Data_QC(:,VAR_N) = SoilMC;
else % There is no data for this variable
    Data_QC(:,VAR_N) = NaN; 
end      
clc;

%% MVF - 59 - Soil Moisture D 
disp('Soil Moisture D ')
VAR_N = 59;
Var_file_save = '\SoilMD.mat';

% Check we have any data
if ~any(noData == VAR_N) % There is some data for this variable
    if exist(strcat(QCdir,Var_file_save)) == 2
        disp('Previous Quality Controled Soil Moisture D   down found')
        u_inp = user_inp(cellstr('Enter 1 - Use this data. 2 - Create new QC data.'), 1, 2);
    else
        u_inp = 2;
    end

    if u_inp == 2
        % Get Raw data
        Var_raw = Data_Raw(:,VAR_N); 
        % QC (Interactive)
        [SoilMD, accepted] = data_QC_NW_edits(Var_raw, 1, TIME_out, 10);
        if accepted == 0
            disp('Exiting because QC was aborted')
            return
        end
        % Save Var QC data
        eval(['save ' strcat(QCdir,Var_file_save) ' SoilMD']);
    elseif u_inp == 1
       load(strcat(QCdir,Var_file_save)) 
    end
    Data_QC(:,VAR_N) = SoilMD;
else % There is no data for this variable
    Data_QC(:,VAR_N) = NaN; 
end      
clc;

%% MVF - 60 - Soil Moisture E 
disp('Soil Moisture E ')
VAR_N = 60;
Var_file_save = '\SoilME.mat';

% Check we have any data
if ~any(noData == VAR_N) % There is some data for this variable
    if exist(strcat(QCdir,Var_file_save)) == 2
        disp('Previous Quality Controled Soil Moisture E   down found')
        u_inp = user_inp(cellstr('Enter 1 - Use this data. 2 - Create new QC data.'), 1, 2);
    else
        u_inp = 2;
    end

    if u_inp == 2
        % Get Raw data
        Var_raw = Data_Raw(:,VAR_N); 
        % QC (Interactive)
        [SoilME, accepted] = data_QC_NW_edits(Var_raw, 1, TIME_out, 10);
        if accepted == 0
            disp('Exiting because QC was aborted')
            return
        end
        % Save Var QC data
        eval(['save ' strcat(QCdir,Var_file_save) ' SoilME']);
    elseif u_inp == 1
       load(strcat(QCdir,Var_file_save)) 
    end
    Data_QC(:,VAR_N) = SoilME;
else % There is no data for this variable
    Data_QC(:,VAR_N) = NaN; 
end      
clc;

%% MVF - 52 - Snow Layer Temperature A
disp('Snow Layer Temperature A ')
VAR_N = 52;
Var_file_save = '\QC_SLT_A.mat';

% Check we have any data
if ~any(noData == VAR_N) % There is some data for this variable
    if exist(strcat(QCdir,Var_file_save)) == 2
        disp('Previous Quality Controled Snow Layer Temperature A  down found')
        u_inp = user_inp(cellstr('Enter 1 - Use this data. 2 - Create new QC data.'), 1, 2);
    else
        u_inp = 2;
    end

    if u_inp == 2
        % Get Raw data
        Var_raw = Data_Raw(:,VAR_N); 
        % QC (Interactive)
        [QC_SLT_A, accepted] = data_QC_NW_edits(Var_raw, 1, TIME_out, 1);
        if accepted == 0
            disp('Exiting because QC was aborted')
            return
        end
        % Save Var QC data
        eval(['save ' strcat(QCdir,Var_file_save) ' QC_SLT_A']);
    elseif u_inp == 1
       load(strcat(QCdir,Var_file_save)) 
    end
    Data_QC(:,VAR_N) = QC_SLT_A;
else % There is no data for this variable
    Data_QC(:,VAR_N) = NaN; 
end      
clc;

%% MVF - 23 - Lysimeter outflow
disp('Lysimeter')
VAR_N = 23;
Var_file_save = '\QC_Lys.mat';

% Check we have any data
if ~any(noData == VAR_N) % There is some data for this variable
    if exist(strcat(QCdir,Var_file_save)) == 2
        disp('Previous Quality Controled Lysimeter data found')
        u_inp = user_inp(cellstr('Enter 1 - Use this data. 2 - Create new QC data.'), 1, 2);
    else
        u_inp = 2;
    end

    if u_inp == 2
        % Get Raw data
        Var_raw = Data_Raw(:,VAR_N); 
        % data_QC wants mm not m
        Var_raw = Var_raw .* 1000; 
        % QC (Interactive)
        [QC_Lys, accepted] = data_QC_NW_edits(Var_raw, 1, TIME_out, 2);
        if accepted == 0
            disp('Exiting because QC was aborted')
            return
        end
        % Convert back to m
        QC_Lys = QC_Lys ./ 1000; 
        % Save Var QC data
        eval(['save ' strcat(QCdir,Var_file_save) ' QC_Lys']);
    elseif u_inp == 1
       load(strcat(QCdir,Var_file_save)) 
    end
    Data_QC(:,VAR_N) = QC_Lys;
else % There is no data for this variable
    Data_QC(:,VAR_N) = NaN; 
end      
clc;

%% MVF - 18 - Leaf Wetness Sensor mV
disp('Leaf Wetness Sensor mV')
VAR_N = 18;
Var_file_save = '\QC_LWS.mat';

% Check we have any data
if ~any(noData == VAR_N) % There is some data for this variable
    if exist(strcat(QCdir,Var_file_save)) == 2
        disp('Previous Quality Controled Leaf Wetness Sensor Data found')
        u_inp = user_inp(cellstr('Enter 1 - Use this data. 2 - Create new QC data.'), 1, 2);
    else
        u_inp = 2;
    end

    if u_inp == 2
        % Get Raw data
        Var_raw = Data_Raw(:,VAR_N);
        % QC (Interactive)
        [QC_LWS, accepted] = data_QC_NW_edits(Var_raw, 1, TIME_out, 8);
        if accepted == 0
            disp('Exiting because QC was aborted')
            return
        end
        % Save Var QC data
        eval(['save ' strcat(QCdir,Var_file_save) ' QC_LWS']);
    elseif u_inp == 1
       load(strcat(QCdir,Var_file_save)) 
    end
    Data_QC(:,VAR_N) = QC_LWS;
else % There is no data for this variable
    Data_QC(:,VAR_N) = NaN; 
end 
clc;


%% MVF - 20 - Snow Depth m
disp('Snow Depth A m')
VAR_N = 20;
Var_file_save = '\QC_SD.mat';

% Check we have any data
if ~any(noData == VAR_N) % There is some data for this variable
    if exist(strcat(QCdir,Var_file_save)) == 2
        disp('Previous Quality Controled Snow Depth Data A found')
        u_inp = user_inp(cellstr('Enter 1 - Use this data. 2 - Create new QC data.'), 1, 2);
    else
        u_inp = 2;
    end

    if u_inp == 2
        % Get Raw data
        Var_raw = Data_Raw(:,VAR_N).*1000;
        % QC (Interactive)
        [QC_SD, accepted] = data_QC_NW_edits(Var_raw, 1, TIME_out, 6);
        QC_SD = QC_SD ./ 1000; % mm to m
        if accepted == 0
            disp('Exiting because QC was aborted')
            return
        end
        % Save Var QC data
        eval(['save ' strcat(QCdir,Var_file_save) ' QC_SD']);
    elseif u_inp == 1
       load(strcat(QCdir,Var_file_save)) 
    end
    Data_QC(:,VAR_N) = QC_SD;
else % There is no data for this variable
    Data_QC(:,VAR_N) = NaN; 
end 
clc;

%% MVF - 24 - Snow Depth  B m
disp('Snow Depth B m')
VAR_N = 24;
Var_file_save = '\QC_SD_B.mat';

% Check we have any data
if ~any(noData == VAR_N) % There is some data for this variable
    if exist(strcat(QCdir,Var_file_save)) == 2
        disp('Previous Quality Controled Snow Depth B Data found')
        u_inp = user_inp(cellstr('Enter 1 - Use this data. 2 - Create new QC data.'), 1, 2);
    else
        u_inp = 2;
    end

    if u_inp == 2
        % Get Raw data
        Var_raw = Data_Raw(:,VAR_N).*1000;
        % QC (Interactive)
        [QC_SD_B, accepted] = data_QC_NW_edits(Var_raw, 1, TIME_out, 6);
        QC_SD_B = QC_SD_B ./ 1000; % mm to m
        if accepted == 0
            disp('Exiting because QC was aborted')
            return
        end
        % Save Var QC data
        eval(['save ' strcat(QCdir,Var_file_save) ' QC_SD_B']);
    elseif u_inp == 1
       load(strcat(QCdir,Var_file_save)) 
    end
    Data_QC(:,VAR_N) = QC_SD_B;
else % There is no data for this variable
    Data_QC(:,VAR_N) = NaN; 
end 
clc;



%% MVF - 48 - Snow Depth 24 change 6am
disp('Snow Depth 24 change from 6am m')
VAR_N = 48;
Var_file_save = '\QC_SD_24hr_6am.mat';

% Check we have any data
if ~any(noData == VAR_N) % There is some data for this variable
    if exist(strcat(QCdir,Var_file_save)) == 2
        disp('Previous Quality Controled Snow Depth 24 change from 6am Data found')
        u_inp = user_inp(cellstr('Enter 1 - Use this data. 2 - Create new QC data.'), 1, 2);
    else
        u_inp = 2;
    end

    if u_inp == 2
        % Get Raw data
        Var_raw = Data_Raw(:,VAR_N).*1000;
        % QC (Interactive)
        [QC_SD_24hr_6am, accepted] = data_QC_NW_edits(Var_raw, 1, TIME_out, 6);
        QC_SD_24hr_6am = QC_SD_24hr_6am ./ 1000; % mm to m
        if accepted == 0
            disp('Exiting because QC was aborted')
            return
        end
        % Save Var QC data
        eval(['save ' strcat(QCdir,Var_file_save) ' QC_SD_24hr_6am']);
    elseif u_inp == 1
       load(strcat(QCdir,Var_file_save)) 
    end
    Data_QC(:,VAR_N) = QC_SD_24hr_6am;
else % There is no data for this variable
    Data_QC(:,VAR_N) = NaN; 
end 
clc;

%% MVF - 49 - 24hr SWE change 6am
disp('24hr SWE change from 6am')
VAR_N = 49;
Var_file_save = '\QC_SWE_24hr_6am.mat';

% Check we have any data
if ~any(noData == VAR_N) % There is some data for this variable
    if exist(strcat(QCdir,Var_file_save)) == 2
        disp('Previous Quality Controled 24hr SWE change from 6am Data found')
        u_inp = user_inp(cellstr('Enter 1 - Use this data. 2 - Create new QC data.'), 1, 2);
    else
        u_inp = 2;
    end

    if u_inp == 2
        % Get Raw data
        Var_raw = Data_Raw(:,VAR_N).*1000;
        % QC (Interactive)
        [QC_SWE_24hr_6am, accepted] = data_QC_NW_edits(Var_raw, 1, TIME_out, 6);
        QC_SWE_24hr_6am = QC_SWE_24hr_6am ./ 1000; % mm to m
        if accepted == 0
            disp('Exiting because QC was aborted')
            return
        end
        % Save Var QC data
        eval(['save ' strcat(QCdir,Var_file_save) ' QC_SWE_24hr_6am']);
    elseif u_inp == 1
       load(strcat(QCdir,Var_file_save)) 
    end
    Data_QC(:,VAR_N) = QC_SWE_24hr_6am;
else % There is no data for this variable
    Data_QC(:,VAR_N) = NaN; 
end 
clc;


%% MVF - 50 - 24hr Density 6am
disp('24hr Density change from 6am')
VAR_N = 50;
Var_file_save = '\QC_DEN_24hr_6am.mat';

% Check we have any data
if ~any(noData == VAR_N) % There is some data for this variable
    if exist(strcat(QCdir,Var_file_save)) == 2
        disp('Previous Quality Controled 24hr SWE change from 6am Data found')
        u_inp = user_inp(cellstr('Enter 1 - Use this data. 2 - Create new QC data.'), 1, 2);
    else
        u_inp = 2;
    end

    if u_inp == 2
        % Get Raw data
        % QC (Interactive)
        [QC_DEN_24hr_6am, accepted] = data_QC_NW_edits(Var_raw, 1, TIME_out, 3);
        if accepted == 0
            disp('Exiting because QC was aborted')
            return
        end
        % Save Var QC data
        eval(['save ' strcat(QCdir,Var_file_save) ' QC_DEN_24hr_6am']);
    elseif u_inp == 1
       load(strcat(QCdir,Var_file_save)) 
    end
    Data_QC(:,VAR_N) = QC_DEN_24hr_6am;
else % There is no data for this variable
    Data_QC(:,VAR_N) = NaN; 
end 
clc;


%% MVF - 44 - Height of CSAT m
disp('Height of CSAT m')
VAR_N = 44;
Var_file_save = '\QC_HSCAT.mat';

% Check we have any data
if ~any(noData == VAR_N) % There is some data for this variable
    if exist(strcat(QCdir,Var_file_save)) == 2
        disp('Previous Quality Controled Height of CSAT Data found')
        u_inp = user_inp(cellstr('Enter 1 - Use this data. 2 - Create new QC data.'), 1, 2);
    else
        u_inp = 2;
    end

    if u_inp == 2
        % Get Raw data
        Var_raw = Data_Raw(:,VAR_N);
        % QC (Interactive)
        [QC_HSCAT, accepted] = data_QC_NW_edits(Var_raw, 1, TIME_out, 6);
        if accepted == 0
            disp('Exiting because QC was aborted')
            return
        end
        % Save Var QC data
        eval(['save ' strcat(QCdir,Var_file_save) ' QC_HSCAT']);
    elseif u_inp == 1
       load(strcat(QCdir,Var_file_save)) 
    end
    Data_QC(:,VAR_N) = QC_HSCAT;
else % There is no data for this variable
    Data_QC(:,VAR_N) = NaN; 
end 
clc;



%% MVF - 45 - CSAT Azimuth from true north
disp('CSAT Azimuth from true north')
VAR_N = 45;
Var_file_save = '\QC_AzCSAT.mat';

% Check we have any data
if ~any(noData == VAR_N) % There is some data for this variable
    if exist(strcat(QCdir,Var_file_save)) == 2
        disp('Previous Quality Controled Height of CSAT Data found')
        u_inp = user_inp(cellstr('Enter 1 - Use this data. 2 - Create new QC data.'), 1, 2);
    else
        u_inp = 2;
    end

    if u_inp == 2
        % Get Raw data
        Var_raw = Data_Raw(:,VAR_N);
        % QC (Interactive)
        [QC_AzCSAT, accepted] = data_QC_NW_edits(Var_raw, 1, TIME_out, 6);
        if accepted == 0
            disp('Exiting because QC was aborted')
            return
        end
        % Save Var QC data
        eval(['save ' strcat(QCdir,Var_file_save) ' QC_AzCSAT']);
    elseif u_inp == 1
       load(strcat(QCdir,Var_file_save)) 
    end
    Data_QC(:,VAR_N) = QC_AzCSAT;
else % There is no data for this variable
    Data_QC(:,VAR_N) = NaN; 
end 
clc;


%% MVF - 47 - Snow Water Equivelent (SWE)
disp('Snow Water Equivelent (SWE)')
VAR_N = 47;
Var_file_save = '\QC_SWE.mat';

% Check we have any data
if ~any(noData == VAR_N) % There is some data for this variable
    if exist(strcat(QCdir,Var_file_save)) == 2
        disp('Previous Quality Controled Snow Water Equivelent Data found')
        u_inp = user_inp(cellstr('Enter 1 - Use this data. 2 - Create new QC data.'), 1, 2);
    else
        u_inp = 2;
    end

    if u_inp == 2
        % Get Raw data
        Var_raw = Data_Raw(:,VAR_N);
        Var_raw = Var_raw.*1000; % m to mm
        % QC (Interactive)
        [QC_SWE, accepted] = data_QC_NW_edits(Var_raw, 1, TIME_out, 6);
        QC_SWE = QC_SWE./1000; % mm to m
        if accepted == 0
            disp('Exiting because QC was aborted')
            return
        end
        % Save Var QC data
        eval(['save ' strcat(QCdir,Var_file_save) ' QC_SWE']);
    elseif u_inp == 1
       load(strcat(QCdir,Var_file_save)) 
    end
    Data_QC(:,VAR_N) = QC_SWE;
else % There is no data for this variable
    Data_QC(:,VAR_N) = NaN; 
end 
clc;

%% MVF - 53 - Freazing level A
disp('Freazing level A')
VAR_N = 53;
Var_file_save = '\QC_FL.mat';

% Check we have any data
if ~any(noData == VAR_N) % There is some data for this variable
    if exist(strcat(QCdir,Var_file_save)) == 2
        disp('Previous Quality Controled Freazing level A Data found')
        u_inp = user_inp(cellstr('Enter 1 - Use this data. 2 - Create new QC data.'), 1, 2);
    else
        u_inp = 2;
    end

    if u_inp == 2
        % Get Raw data
        Var_raw = Data_Raw(:,VAR_N);
        % QC (Interactive)
        [QC_FL, accepted] = data_QC_NW_edits(Var_raw, 1, TIME_out, 6);
        if accepted == 0
            disp('Exiting because QC was aborted')
            return
        end
        % Save Var QC data
        eval(['save ' strcat(QCdir,Var_file_save) ' QC_FL']);
    elseif u_inp == 1
       load(strcat(QCdir,Var_file_save)) 
    end
    Data_QC(:,VAR_N) = QC_FL;
else % There is no data for this variable
    Data_QC(:,VAR_N) = NaN; 
end 
clc;


%% MVF - 54 - RamPen Test A
disp('RamPen Test A')
VAR_N = 54;
Var_file_save = '\QC_RPT.mat';

% Check we have any data
if ~any(noData == VAR_N) % There is some data for this variable
    if exist(strcat(QCdir,Var_file_save)) == 2
        disp('Previous Quality Controled RamPen Test A Data found')
        u_inp = user_inp(cellstr('Enter 1 - Use this data. 2 - Create new QC data.'), 1, 2);
    else
        u_inp = 2;
    end

    if u_inp == 2
        % Get Raw data
        Var_raw = Data_Raw(:,VAR_N);
        Var_raw = Var_raw.*1000; % m to mm
        % QC (Interactive)
        [QC_RPT, accepted] = data_QC_NW_edits(Var_raw, 1, TIME_out, 6);
        QC_RPT = QC_RPT./1000; % mm to m
        if accepted == 0
            disp('Exiting because QC was aborted')
            return
        end
        % Save Var QC data
        eval(['save ' strcat(QCdir,Var_file_save) ' QC_RPT']);
    elseif u_inp == 1
       load(strcat(QCdir,Var_file_save)) 
    end
    Data_QC(:,VAR_N) = QC_RPT;
else % There is no data for this variable
    Data_QC(:,VAR_N) = NaN; 
end 
clc;


%% Create matrix of missing timesteps (for validating model results, etc.)
Missing_I = isnan(Data_QC); 

%% Save QC data
disp('Saving QC Data')
file_name_out = strcat(QCdir,'\',file_out);
eval(['save ' file_name_out ' Data_QC TIME_out Missing_I MHeight_MVF site_E'])
