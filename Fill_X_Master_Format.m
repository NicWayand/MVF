clear all; close all; clc; smp;
tic

% This script takes QC data in MVF (Created by
% QC_X_Master_Formate.m), QC_RHs according to user specifications, 
% saves a continous .mat file in MVF
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

%% User Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
site = 'SNQ15';

%% Input info
Forcing_N_in  = 'OctMar';
file_in = strcat(site,'_QC.mat');
dt_in = 0.5; 

%% Output info
Forcing_N_out = 'OctMar';
file_out = strcat(site,'_Filled.mat');

% Spinup Settings (option to add spin up period to forcing file)
Add_Spinup = false;

% Other info for cloud factor estimationg (only needed if DOKIA LW
% esitmation is used

% MFS
%latitude  = 39.94855;
%longitude = -105.196061; 
%UTCoffset = -7; % only used for LW estimation

% SNQ
% latitude  = 47.424865;
% longitude = -121.413831; 
% UTCoffset = -8; % only used for LW estimation

% Name of NLDAS FILE (comment out if not used)
% LDAS_file = 'SNQ_WY2013'; % Wy 2013 Complete
% LDAS_file = 'SNQ_2014_Oct_May.mat';
% LDAS_file = ''; % Oct 2012 to Present


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Directories
homedir = 'G:\ncar\d1\data\Point_Snow_Sites\'; % REPLACE YOUR HOME DIR HERE
QCdir = fullfile(homedir,'QC',site,Forcing_N_in);
Filldir = fullfile(homedir,'Fill',site,Forcing_N_out);
if exist(Filldir,'dir') ~= 7
    mkdir(Filldir)
end
cd(Filldir)

%% Load NLDAS data
if exist('LDAS_file','var')
    load(['G:\ncar\d1\data\LDAS\matlabF\' LDAS_file]) % REPLACE YOUR HOME DIR HERE
end
%% Load QC MVF data
load(strcat(QCdir,'\',file_in)); 
TIME_out_all = time_builder(TIME_out);

%% Filled Input Data
Data_Fill = NaN(size(Data_QC)); % QC'ed data (only includes variables QC'ed below)

%% MVF - 1 - Air Temperature A (assumed the best temperature data here, options allow use of other temperature data to fill)
disp('Air Temperature')
VAR_N = 1; 
alternate_VAR_N = [31 42 43];
PriorityID = ['A';'B';'C';'D'];
Var_file_save = '\Fill_Temp.mat';

if exist(strcat(Filldir,Var_file_save)) == 2
    disp('Previous Quality Controled Temp found')
    u_inp = user_inp(cellstr('Enter 1 - Use this data. 2 - Create new Fill data. 3 - Continue with previous filled data'), 1, 2, 3);
else
    u_inp = 2;
end

if u_inp == 2 || u_inp == 3
    % check for other Temperature data (i.e. B,C,D)
    I_alt_data = sum(~isnan(Data_QC(:,alternate_VAR_N)),1) > 0;
    if any(I_alt_data) % Found some other data
        % grab other data
        fprintf('Found %d alternate data sets\n',sum(I_alt_data))
        Var_raw_alt = Data_QC(:,alternate_VAR_N);
    end
    
    % Check for NLDAS data
    if exist('LDAS_file','var')
        fprintf('Found NLDAS temperature data, use? (will neglect other alternates)\n')
        u_inpnldas = user_inp(cellstr('Enter 1 - Use NLDAS. 2 - Do not use'), 1, 2);
    end
 
    if u_inp == 2
        % Get QC data
        Var_raw = Data_QC(:,VAR_N);
    elseif u_inp == 3 % Get previous Filled data
        load(strcat(Filldir,Var_file_save)) 
        Var_raw = Fill_Temp;    
    end
    
    if u_inpnldas == 1
        % Find similar LDAS SW  and Data_Fill times
        [C,IA,IB] = intersect(datenum_round_off(TIME_out,'minute'),datenum_round_off(TIME_30min(:,7),'minute'));

        Temp_LDAS = data_NLDAS(IB,3); 
        Temp_LDAS = Temp_LDAS - 273.15; % K to C
        
        [Fill_Temp_all] = data_fill_NW_edits(TIME_out_all(IA,:), [Var_raw(IA) Temp_LDAS], 1, 40, -30, 1, 1, 48, [site_E site_E+1]);  
        Fill_Temp(IA) = Fill_Temp_all(:,1); % Grab first filled data set (original 1)
        
    elseif any(I_alt_data) % Include other data sets for filling
                
        disp('Now plotting all Temperature data found')
        figure(1)
        hold on
        p1 = plot(TIME_out_all(:,7),Var_raw,'k','linewidth',3);
        p2 = plot(TIME_out_all(:,7),Var_raw_alt);
        tlabel
        legend([p1 p2(:)'],PriorityID)
        pause
        
        % Fill (Interactive)
        [Fill_Temp_all] = data_fill_NW_edits(TIME_out_all, [Var_raw Var_raw_alt], 1, 40, -30, 1, 1, 48, site_E + (1:size(Var_raw_alt,2)+1) );  
        Fill_Temp = Fill_Temp_all(:,1); % Grab first filled data set (original 1)
    else % Dont use other data 
        % Fill (Interactive) (one data set)
        [Fill_Temp] = data_fill_NW_edits(TIME_out_all, Var_raw, 1, 40, -30, 1, 1, 48, site_E);
    end
    % Save Var Fill data
    eval(['save ' strcat(Filldir,Var_file_save) ' Fill_Temp']);
elseif u_inp == 1
   load(strcat(Filldir,Var_file_save)) 
end
Data_Fill(:,VAR_N) = Fill_Temp;
clc; clear alternate_VAR_N

%% MVF - 2 - Precipitation
disp('Precipitation')
VAR_N = 2;
alternate_VAR_N = [49 55];
Var_file_save = '\Fill_Precip.mat';

if exist(strcat(Filldir,Var_file_save)) == 2
    disp('Previous Quality Controled Precip found')
    u_inp = user_inp(cellstr('Enter 1 - Use this data. 2 - Create new Fill data. 3 - Continue with previous filled data'), 1, 2, 3);
else
    u_inp = 2;
end

if u_inp == 2 || u_inp == 3
    % check for other Temperature data (i.e. B,C,D)
    I_alt_data = sum(~isnan(Data_QC(:,alternate_VAR_N)),1) > 0;
    if any(I_alt_data) % Found some other data
        % grab other data
        fprintf('Found %d alternate data sets\n',sum(I_alt_data))
        Var_raw_alt = Data_QC(:,alternate_VAR_N);
        Var_raw_alt = Var_raw_alt .* 1000;
    end
    
    % Check for NLDAS data
    if exist('LDAS_file','var')
        fprintf('Found NLDAS temperature data, use? (will neglect other alternates)\n')
        u_inpnldas = user_inp(cellstr('Enter 1 - Use NLDAS. 2 - Do not use'), 1, 2);
    end
    
    if u_inp == 2
        % Get QC data
        Var_raw = Data_QC(:,VAR_N);
    elseif u_inp == 3 % Get previous Filled data
        load(strcat(Filldir,Var_file_save)) 
        Var_raw = Fill_Precip;    % data_fill wants mm not m
    end
    % Convert from m to mm
    Var_raw = Var_raw .* 1000; 
    
    if u_inpnldas == 1
        % Find similar LDAS SW and Data_Fill times
        [C,IA,IB] = intersect(datenum_round_off(TIME_out,'minute'),datenum_round_off(TIME_30min(:,7),'minute')); 

        Precip_LDAS = data_NLDAS(IB,5); 
        Precip_LDAS = kg_per_m2_to_m(Precip_LDAS,((data_NLDAS(IB,3)<0)+1) ); % kg/m^3 to m/hr.... ((data_NLDAS(IB,3)<0)+1) defines what precip type (roughly)
        Precip_LDAS = Precip_LDAS .*1000; % m to mm
        
        
        figure(1)
        hold on
        plot(C,(Var_raw(IA)),'b')
        plot(C,(Precip_LDAS),'k')
        ylabel('Incremental precip mm')
        legend('Observed','NLDAS')
        tlabel
        pause
        close(1)
        
        figure(2)
        hold on
        plot(C,nancumsum(Var_raw(IA)),'b')
        plot(C,nancumsum(Precip_LDAS),'k')
        ylabel('Cummulative precip mm')
        legend('Observed','NLDAS')
        tlabel
        pause
        close(2)
        
        % Fill (Interactive) with NLDAS data
        [Fill_Precip_all] = data_fill_NW_edits(TIME_out_all(IA,:), [Var_raw(IA) Precip_LDAS], 2, 100, 0, 1, 1, 48, [site_E site_E+1]);
        Fill_Precip_all = Fill_Precip_all ./ 1000; % Convert back from mm to m 
        Fill_Precip       = Var_raw./ 1000; % first initilize to orig values, mm to m 
        Fill_Precip(IA)   = Fill_Precip_all(:,1); % Grab first filled data set (original 1), assinged new update values
        
    elseif any(I_alt_data) % Include other data sets for filling
        
        % either have 24 SWE or other incremental precip data
        if alternate_VAR_N(I_alt_data) == 49 % 24 SWE
            fprintf('Found 24 hours SWE values, Making sure Precip is >= 24 SWE accumlation\n')
            Var_raw_alt = Var_raw_alt(:,I_alt_data);
            
            % Agg to Daily values
            fprintf('Assuming SWE value is accumlation at 6 am\n')
            Ifirst6am = find(TIME_out_all(:,4) == 6 & TIME_out_all(:,5) == 0,1,'first'); % Do not include precip at 6-6:30am. 
            [Var_raw_daily,TIME_temp_dly] = aggMET(Var_raw,TIME_out_all,dt_in,24,TIME_out_all(Ifirst6am,1),TIME_out_all(Ifirst6am,2),TIME_out_all(Ifirst6am,3),TIME_out_all(Ifirst6am,4),0,1,5,100);
            Ifirst6am_24hrs = find(TIME_out_all(:,4) == 6 & TIME_out_all(:,5) == 0,1,'first') + 2; % Do not include precip at 6-6:30am. 
            [Var_raw_alt_daily,~] = aggMET(Var_raw_alt,TIME_out_all,dt_in,24,TIME_out_all(Ifirst6am_24hrs,1),TIME_out_all(Ifirst6am_24hrs,2),TIME_out_all(Ifirst6am_24hrs,3),TIME_out_all(Ifirst6am_24hrs,4),0,1,5,100);

            % Set observed precip data days with all NaN, to zero so that
            % difference from daily SWE accums will be positive and can be
            % added to precip (distributed is some way)
            Var_raw_daily(isnan(Var_raw_daily)) = 0;

            % Calculate difference
            Daily_Accum_Diff = Var_raw_daily - Var_raw_alt_daily; 

            figure; hold on
            plot(TIME_temp_dly(:,7),Var_raw_daily,'k')
            plot(TIME_temp_dly(:,7),Var_raw_alt_daily,'r')

            figure; hold on
            plot(TIME_temp_dly(:,7),Daily_Accum_Diff)
            legend('24hrs Obs Precip - 24hrs Obs SWE')
            tlabel

            pause

            % Find days when Precip is less then SWE
            I_bad = Daily_Accum_Diff < 0;

            fprintf('Found %f mm missed by the Precip Guage in total\n',sum(Daily_Accum_Diff(I_bad)))

            fprintf('Do you want to bias correct the Precip Gauge?\n')
            apply_Precip_correction = user_inp(cellstr('Enter 1 - Yes. 2 - No'), 1, 2);
            if apply_Precip_correction
                % Set good ratio days to zero
                Daily_Accum_Diff(~I_bad) = 0; 

                % Make negative differences positive
                Daily_Accum_Diff = Daily_Accum_Diff * -1;

                fprintf('How do you want to disaggragate 24 hour biases to fine timestep?\n')
                OPTION_Precip_correction = user_inp(cellstr('Enter 1 - Uniformly. 3 - Timesteps with precipitation already'), 1, 3);

                % Disaggragate new Precip to original time step
                Var_raw_corr = disaggMET(TIME_temp_dly(:,7),Daily_Accum_Diff,24,dt_in,TIME_out,1,OPTION_Precip_correction,0,Var_raw);



                figure; hold on
                plot(TIME_out,Var_raw,'k')

                % Add to origial value 
                Var_new = Var_raw + Var_raw_corr;
                % Add corrected precip when orig is Nan and Var_raw_corr has
                % data
                I_N = isnan(Var_raw) & Var_raw_corr ~= 0;
                Var_new(I_N) = Var_raw_corr(I_N);

                plot(TIME_out,Var_new,'r')
                legend('Orig','With correction')
                tlabel

                if nansum(Var_new) - nansum(nansum([Var_raw Var_raw_corr])) > 10^06
                    error('missing some precip')
                end

                pause


                fprintf('Do to use the corrected Preecip data?\n')
                use_corr_Precip = user_inp(cellstr('Enter 1 - Yes. 2 - No'), 1, 2);
                if use_corr_Precip
                    Fill_Precip = Var_new;
                else
                    Fill_Precip = Var_raw;
                end

                Fill_Precip = Fill_Precip ./ 1000; % mm to m

            else
                disp('Alt Precip not used')
            end
        elseif alternate_VAR_N(I_alt_data) == 55 % Alternate incremental precip data
            Var_raw_alt = Var_raw_alt(:,I_alt_data);

            % Plot two data sets
            figure(1)
            hold on
            plot(TIME_out,Var_raw,'b')
            plot(TIME_out,Var_raw_alt,'k')
            ylabel('Incremental precip mm')
            legend('Precip A','Precip B')
            tlabel
            pause
            close(1)
            
%             % Calc total bias for alternate data source of precip using available
%             % precip MVF = 2
%             disp('Calc. Bias between sources')
%             disp('WARNING! Asssuming same time step (i.e. both daily)')
%             Var_alt_bias = nanmean(Var_raw_alt-Var_raw);
%             Var_raw_alt = Var_raw_alt + Var_alt_bias;
%             fprintf('Bias is %f\n\n',Var_alt_bias)
%             
%             disp('Now plotting with bias corrected')
%             % Plot two data sets
%             figure(1)
%             hold on
%             plot(TIME_out,Var_raw,'b')
%             plot(TIME_out,Var_raw_alt,'k')
%             ylabel('Incremental precip mm')
%             legend('Precip A','Precip B')
%             tlabel
%             pause
%             close(1)

            
            
            % Fill (Interactive) with NLDAS data
            [Fill_Precip_all] = data_fill_NW_edits(TIME_out_all, [Var_raw Var_raw_alt], 2, 100, 0, 1, 1, 48, [site_E site_E+1]);
            Fill_Precip  = Fill_Precip_all(:,1) ./ 1000; % Convert back from mm to m 

        end
        
    else
        % Fill (Interactive)
        [Fill_Precip] = data_fill_NW_edits(TIME_out_all, Var_raw, 2, 500, 0, 1, 1, 48, site_E);
        Fill_Precip = Fill_Precip ./ 1000; % mm to m
    end
    
    % Save Var Fill data
    eval(['save ' strcat(Filldir,Var_file_save) ' Fill_Precip']);
elseif u_inp == 1
   load(strcat(Filldir,Var_file_save)) 
end
Data_Fill(:,VAR_N) = Fill_Precip;
clc; clear alternate_VAR_N

%% MVF - 4 - Relative Humidity
disp('Relative Humidity')
VAR_N = 4;
alternate_VAR_N = [32];
PriorityID = ['A';'B';'C';'D'];
Var_file_save = '\Fill_RH.mat';

if exist(strcat(Filldir,Var_file_save)) == 2
    disp('Previous Quality Controled Relative Humidity found')
    u_inp = user_inp(cellstr('Enter 1 - Use this data. 2 - Create new Fill data. 3 - Continue with previous filled data'), 1, 2, 3);
else
    u_inp = 2;
end

if u_inp == 2 || u_inp == 3
     % check for other RH data (i.e. B,C,D)
    I_alt_data = sum(~isnan(Data_QC(:,alternate_VAR_N)),1) > 0;
    if any(I_alt_data) % Found some other data
        % grab other data
        fprintf('Found %d alternate data sets\n',sum(I_alt_data))
        Var_raw_alt = Data_QC(:,alternate_VAR_N);
        Var_raw_alt = Var_raw_alt ./ 100; 

    end
    
    % Check for NLDAS data
    if exist('LDAS_file','var')
        fprintf('Found NLDAS temperature data, use? (will neglect other alternates)\n')
        u_inpnldas = user_inp(cellstr('Enter 1 - Use NLDAS. 2 - Do not use'), 1, 2);
    end
    
    if u_inp == 2
        % Get raw QC data
        Var_raw = Data_QC(:,VAR_N);
    elseif u_inp == 3 % Get previous Filled data
        load(strcat(Filldir,Var_file_save)) 
        Var_raw = Fill_RH;
    end
    % Convert from % to fractional
    Var_raw     = Var_raw ./ 100; 
    
    if u_inpnldas == 1
        % Find similar LDAS SW and Data_Fill times
        [C,IA,IB] = intersect(datenum_round_off(TIME_out,'minute'),datenum_round_off(TIME_30min(:,7),'minute')); 

        SH_LDAS = data_NLDAS(IB,4); 
        RH_LDAS = Q_to_RH(SH_LDAS,data_NLDAS(IB,3),hPa_to_Pa(data_NLDAS(IB,7)));
        RH_LDAS = RH_LDAS ./ 100; 
        
        figure(1)
        hold on
        plot(C,Var_raw(IA),'k','linewidth',3)
        plot(C,RH_LDAS)
        tlabel
        legend('Obs','NLDAS')
        pause
        close(1)
        
        [Fill_RH_all] = data_fill_NW_edits(TIME_out_all(IA,:), [Var_raw(IA) RH_LDAS], 3, 1, 0, 1, 1, 48, [site_E site_E+1]);  
        Fill_RH(IA) = Fill_RH_all(:,1); 
    
    elseif any(I_alt_data) % Include other data sets for filling
        disp('Now plotting all RH data found')
        figure(1)
        hold on
        p1 = plot(TIME_out_all(:,7),Var_raw,'k','linewidth',3);
        p2 = plot(TIME_out_all(:,7),Var_raw_alt);
        tlabel
        legend([p1 p2(:)'],PriorityID)
        pause
        close(1)

        % Fill (Interactive)
        [Fill_RH_all] = data_fill_NW_edits(TIME_out_all, [Var_raw Var_raw_alt], 3, 1, 0, 1, 1, 48, [site_E site_E.*ones(1,size(Var_raw_alt,2))+1 ]);  
        Fill_RH = Fill_RH_all(:,1); 
    else
        % Fill (Interactive) (no other alt data)
        [Fill_RH] = data_fill_NW_edits(TIME_out_all, Var_raw, 3, 1, 0, 1, 1, 48, site_E);
    end

    Fill_RH = Fill_RH .* 100; % Fraction to %
    
    % Double check to remove values greater than 100 (after a lapse rate)
    disp('Forcing all RH output to be within 0-100 range')
    Fill_RH(Fill_RH>100) = 100;
    Fill_RH(Fill_RH<0) = 0; 
    
    % Save Var Fill data
    eval(['save ' strcat(Filldir,Var_file_save) ' Fill_RH']);
elseif u_inp == 1
   load(strcat(Filldir,Var_file_save)) 
end
Data_Fill(:,VAR_N) = Fill_RH;
clc; clear alternate_VAR_N

%% MVF - 5 - Wind Speed A
disp('Wind Speed A')
VAR_N = 5;
alternate_VAR_N = [6];
PriorityID = ['A';'B';'C';'D'];
Var_file_save = '\Fill_WS2m.mat';

if exist(strcat(Filldir,Var_file_save)) == 2
    disp('Previous Filled Wind Speed data found')
    u_inp = user_inp(cellstr('Enter 1 - Use this data. 2 - Create new Fill data. 3 - Continue with previous filled data'), 1, 2, 3);
else
    disp('No previous data found, starting filling procedure')
    u_inp = 2;
end

if u_inp == 2 || u_inp == 3
    % check for other WS data (i.e. B,C,D)
    I_alt_data = sum(~isnan(Data_QC(:,alternate_VAR_N)),1) > 0;
    if any(I_alt_data) % Found some other data
        % grab other data
        fprintf('Found %d alternate data sets\n',sum(I_alt_data))
        Var_raw_alt = Data_QC(:,alternate_VAR_N);
        for caltdata = 1:length(alternate_VAR_N)
            fprintf('Scaling alternative wind speed data from %f to %f\n',MHeight_MVF(alternate_VAR_N(caltdata)),MHeight_MVF(VAR_N))
            Var_raw_alt(:,caltdata) = scale_WS(Var_raw_alt(:,caltdata),MHeight_MVF(alternate_VAR_N(caltdata)),MHeight_MVF(VAR_N),1); % logrithmicly
        end
    end
    
    if u_inp == 2
        % Get raw QC data
        Var_raw = Data_QC(:,VAR_N);
    elseif u_inp == 3 % Get previous Filled data
        load(strcat(Filldir,Var_file_save)) 
        Var_raw = Fill_WS2m;
    end
    
    if any(I_alt_data) % Include other data sets for filling
        disp('Now plotting all WS data found')
        figure(1)
        hold on
        p1 = plot(TIME_out_all(:,7),Var_raw,'k','linewidth',3);
        p2 = plot(TIME_out_all(:,7),Var_raw_alt);
        tlabel
        legend([p1 p2(:)'],PriorityID)
        pause
        close(1)
        disp('Scale alternative wind speed by primary wind speed?')
        u_scale_alt = user_inp(cellstr('Enter 1 - Yes. 2 - No.'), 1, 2);
        if u_scale_alt == 1
            % Calculate scaling factor from alternate to primary wind speed
            fctrWS1 = nanmean(Var_raw./Var_raw_alt);
            % Apply factor
            Var_raw_alt = Var_raw_alt .* fctrWS1;
        end

        % Fill (Interactive)
        [Fill_WS2m_all] = data_fill_NW_edits(TIME_out_all, [Var_raw Var_raw_alt], 4, 30, 0, 1, 1, 48, [site_E site_E.*ones(1,size(Var_raw_alt,2))+1 ]);  
        Fill_WS2m = Fill_WS2m_all(:,2); 
    else
        % Fill (Interactive)
        [Fill_WS2m] = data_fill_NW_edits(TIME_out_all, Var_raw, 4, 30, 0, 1, 1, 48, site_E);
    end
    % Save Var Fill data
    eval(['save ' strcat(Filldir,Var_file_save) ' Fill_WS2m']);
elseif u_inp == 1
   load(strcat(Filldir,Var_file_save)) 
end
Data_Fill(:,VAR_N) = Fill_WS2m;
clc; clear alternate_VAR_N


% %% MVF - 14 - Pressure (adjusted)
% disp('Pressure')
% VAR_N = 14;
% Var_file_save = '\Fill_Press.mat';
% 
% if exist(strcat(Filldir,Var_file_save)) == 2
%     disp('Previous Filled Pressure data found')
%     u_inp = user_inp(cellstr('Enter 1 - Use this data. 2 - Create new Fill data. 3 - Continue with previous filled data'), 1, 2, 3);
% else
%     disp('No previous data found, starting filling procedure')
%     u_inp = 2;
% end
% 
% if u_inp == 2 || u_inp == 3
%     if u_inp == 2
%         % Get raw QC data
%         Var_raw = Data_QC(:,VAR_N);
%     elseif u_inp == 3 % Get previous Filled data
%         load(strcat(Filldir,Var_file_save)) 
%         Var_raw = Fill_Press;
%     end
%     % Fill (Interactive)
%     [Fill_Press] = data_fill_NW_edits(TIME_out_all, Var_raw, 6, 1100, 0, 1, 1, 48, site_E);
% %     [Fill_RH] = data_fill_NW_edits(TIME_out_all, Var_raw, 3, 100, 0, 1, site_E);
%     % Save Var Fill data
%     eval(['save ' strcat(Filldir,Var_file_save) ' Fill_Press']);
% elseif u_inp == 1
%    load(strcat(Filldir,Var_file_save)) 
% end
% Data_Fill(:,VAR_N) = Fill_Press;
% clc; clear alternate_VAR_N

%% MVF - 16 - Unadjusted Pressure
disp('Unadjusted Pressure')
VAR_N = 16;
alternate_VAR_N = [46];
PriorityID = ['A';'B';'C';'D'];
Var_file_save = '\Fill_PressUnAdj.mat';

if exist(strcat(Filldir,Var_file_save)) == 2
    disp('Previous Filled Unadjusted Pressure data found')
    u_inp = user_inp(cellstr('Enter 1 - Use this data. 2 - Create new Fill data. 3 - Continue with previous filled data'), 1, 2, 3);
else
    disp('No previous data found, starting filling procedure')
    u_inp = 2;
end

if u_inp == 2 || u_inp == 3
    
    % check for other Pressure data (i.e. B,C,D)
    I_alt_data = sum(~isnan(Data_QC(:,alternate_VAR_N)),1) > 0;
    if any(I_alt_data) % Found some other data
        % grab other data
        fprintf('Found %d alternate data sets\n',sum(I_alt_data))
        Var_raw_alt = Data_QC(:,alternate_VAR_N);
    end
    
    if u_inp == 2
        % Get raw QC data
        Var_raw = Data_QC(:,VAR_N);
    elseif u_inp == 3 % Get previous Filled data
        load(strcat(Filldir,Var_file_save)) 
        Var_raw = Fill_PressUnAdj;
    end
    
    
    % Allow for alternative Pressure source
    Nnan = sum(isnan(Data_QC(:,VAR_N)))/numel(Data_QC(:,VAR_N)) * 100;
    sprintf('There are %3.1f percent missing values',Nnan)
    disp('Do you wish to fill missing values with alternative Pressure source?')
    lat_inp = user_inp(cellstr('Enter 1 - NO. 2 - YES, use alt data source (if available). 3 - YES, use LDAS data'), 1, 2, 3);
    
    % Fill (Interactive)
    if lat_inp == 1 % Simple filling methods
        [Fill_PressUnAdj] = data_fill_NW_edits(TIME_out_all, Var_raw, 6, 1100, 700, 1, 1, 48, site_E);
    elseif lat_inp == 2 % Fill with alt data
        disp('Now plotting all Pressure data found')
        figure(1)
        hold on
        p1 = plot(TIME_out_all(:,7),Var_raw,'k','linewidth',3);
        p2 = plot(TIME_out_all(:,7),Var_raw_alt);
        tlabel
        legend([p1 p2(:)'],PriorityID)
        pause
        
        [Fill_PressUnAdj_all] = data_fill_NW_edits(TIME_out_all, [Var_raw Var_raw_alt], 6, 1100, 700, 1, 1, 48, [site_E site_E.*ones(1,size(Var_raw_alt,2))+1 ]);
        Fill_PressUnAdj = Fill_PressUnAdj_all(:,1); % Grab origial filled data
    elseif lat_inp == 3 % Fill with LDAS data
        % Find similar LDAS SW  and Data_Fill times
        [C,IA,IB] = intersect(datenum_round_off(TIME_out,'minute'),datenum_round_off(TIME_30min(:,7),'minute')); 

        Press_LDAS = data_NLDAS(IB(1):IB(end),7); 
        
        % Fill with LDAS as second station
        [Fill_PressUnAdj_out] = data_fill_NW_edits(TIME_out_all(IA(1):IA(end),:), [Var_raw(IA(1):IA(end)) Press_LDAS], 6, 1100, 700, 1, 1, 48, [site_E site_E+1]);
        % Grab only the first station (combined timeseries)
        Var_raw(IA(1):IA(end)) = Fill_PressUnAdj_out(:,1); 
        Fill_PressUnAdj = Var_raw; % just for consistent output
        
        % Plot to check
        figure(1)
        hold on
        plot(TIME_out,Fill_PressUnAdj,'k')
        plot(TIME_out,Press_LDAS)
        legend('Combined','NLDAS data')
        tlabel
        pause
        close(1)
                
    end
        
    
        
    % Save Var Fill data
    eval(['save ' strcat(Filldir,Var_file_save) ' Fill_PressUnAdj']);
elseif u_inp == 1
   load(strcat(Filldir,Var_file_save)) 
end
Data_Fill(:,VAR_N) = Fill_PressUnAdj;
clc; clear alternate_VAR_N

%% MVF - 9 - Shortwave down
disp('Shortwave down')
VAR_N = 9;
Var_file_save = '\Fill_SWd.mat';

if exist(strcat(Filldir,Var_file_save)) == 2
    disp('Previous Filled Shortwave data found')
    u_inp = user_inp(cellstr('Enter 1 - Use this data. 2 - Create new Fill data. 3 - Continue with previous filled data'), 1, 2, 3);
else
    disp('No previous data found, starting filling procedure')
    u_inp = 2;
end

if u_inp == 2 || u_inp == 3
    if u_inp == 2
        % Get raw QC data
        Var_raw = Data_QC(:,VAR_N);
    elseif u_inp == 3 % Get previous Filled data
        load(strcat(Filldir,Var_file_save)) 
        Var_raw = Fill_SWd;
    end
    
    % Allow for alternative longwave down source
    Nnan = sum(isnan(Data_QC(:,VAR_N)))/numel(Data_QC(:,VAR_N)) * 100;
    sprintf('There are %3.1f percent missing values',Nnan)
    disp('Do you wish to fill missing values with alternative estimated SW down?')
    lat_inp = user_inp(cellstr('Enter 1 - NO. 2 - YES, use DOKIA method. 3 - YES, use LDAS data'), 1, 2, 3);
    
    
    % Fill (Interactive)
    if lat_inp == 1 % Simple filling methods
        [Fill_SWd] = data_fill_NW_edits(TIME_out_all, Var_raw, 1, 1365, 0, 1, 1, 48, site_E);
        
    elseif lat_inp == 2 % Fill with empirical SW method
        error('not implemented')
    elseif lat_inp == 3 % Fill with LDAS data
        % Find similar LDAS SW  and Data_Fill times
        [C,IA,IB] = intersect(datenum_round_off(TIME_out,'minute'),datenum_round_off(TIME_30min(:,7),'minute')); 
%         ISldas = find(TIME_30min(:,7) == TIME_out(1)); if isempty(ISldas); error('ISldas not found, check NLDAS dates'); end;
%         IEldas = find(TIME_30min(:,7) == TIME_out(end)); if isempty(IEldas); error('ISldas not found, check NLDAS dates'); end;
%         data_NLDAS formate: [SW_30min LW_30min T_30min SH_30min Precip_30min W_30min Press_30min];
        SW_LDAS = data_NLDAS(IB(1):IB(end),1); 
        
        % Fill with LDAS as second station
        [Fill_SWd_out] = data_fill_NW_edits(TIME_out_all(IA(1):IA(end),:), [Var_raw(IA(1):IA(end)) SW_LDAS], 1, 1365, 0, 1, 1, 48, [site_E site_E+1]);
        % Grab only the first station (combined timeseries)
        Var_raw(IA(1):IA(end)) = Fill_SWd_out(:,1); 
        Fill_SWd = Var_raw; % just for consistent output
        
        % Plot to check
        figure(1)
        hold on
        plot(TIME_out,Var_raw,'k')
        plot(TIME_out,Fill_SWd)
        tlabel
        pause
        close(1)
                
    end
        
    disp('WARNING, need to add SW to data_fill_NW_edits')
    % Save Var Fill data
    eval(['save ' strcat(Filldir,Var_file_save) ' Fill_SWd']);
elseif u_inp == 1
   load(strcat(Filldir,Var_file_save)) 
end
Data_Fill(:,VAR_N) = Fill_SWd;
clc; clear alternate_VAR_N




%% MVF - 11 - Longwave down
disp('Longwave down')
VAR_N = 11;
Var_file_save = '\Fill_LWd.mat';

if exist(strcat(Filldir,Var_file_save)) == 2
    disp('Previous Filled Longwave data found')
    u_inp = user_inp(cellstr('Enter 1 - Use this data. 2 - Create new Fill data. 3 - Continue with previous filled data'), 1, 2, 3);
else
    disp('No previous data found, starting filling procedure')
    u_inp = 2;

end

if u_inp == 2 || u_inp == 3
    if u_inp == 2
        % Get raw QC data
        Var_raw = Data_QC(:,VAR_N);
    elseif u_inp == 3 % Get previous Filled data
        load(strcat(Filldir,Var_file_save)) 
        Var_raw = Fill_LWd;
    end
    
    % Allow for alternative longwave down source
    Nnan = sum(isnan(Data_QC(:,VAR_N)))/numel(Data_QC(:,VAR_N)) * 100;
    sprintf('There are %3.1f percent missing values',Nnan)
    disp('Do you wish to fill missing values with alternative estimated LW down?')
    lat_inp = user_inp(cellstr('Enter 1 - NO. 2 - YES, use DOKIA method. 3 - YES, use LDAS data'), 1, 2, 3);
    
    % Fill (Interactive)
    if lat_inp == 1 % Simple filling methods
        [Fill_LWd] = data_fill_NW_edits(TIME_out_all, Var_raw, 1, 500, 100, 1, 1, 48, site_E);
        
    elseif lat_inp == 2 % Fill with DOKIA method's LWd
        % Estimate LW down using DOKIA method 
        
        % Estimate cloud factor
        disp('Estimating cloud factor based on SW down from Data_Fill (column 9)')
        [c,~,~,~,~,~] = cloudfactor2_Jessica(TIME_out_all,latitude,longitude,site_E,UTCoffset,Data_Fill(:,9));

        % Estimate LW
        inc_LW        = inc_LW_dokia(TIME_out_all, Data_Fill(:,1), RH_to_fraction(Data_Fill(:,4)), c');

        % Calc and remove LW Bias from observation points
        dokia_rmse = RMSE(Var_raw,inc_LW);
        disp('Removing Bias from DOKIA estimated LW using all observed LW')
%         dokia_bias = nanmean(Var_raw)-nanmean(inc_LW);
        
        % Plot check
        figure(1)
        hold on
%         plot(Var_raw,inc_LW+dokia_rmse,'k*')
%         r1 = refline(1,0);
        plot(TIME_out,Var_raw,'k')
        plot(TIME_out,inc_LW+dokia_rmse,'b')
        legend('Obs-LW-dw','DOKIA-LW-dw')
        tlabel
        pause
                
        % Fill
        [Fill_LWd] = data_fill_NW_edits(TIME_out_all, [Var_raw inc_LW+dokia_rmse], 1, 500, 100, 1, 1, 48, [site_E site_E+1]);
                
        
        disp('Do want to use the filled data set, or ALL DOKIA?')
        lat_inp_2 = user_inp(cellstr('Enter 1 - Filled. 2 - ALL DOKIA'), 1, 2);
        if lat_inp_2 == 1
            % Grab only the first station (combined timeseries)
            Fill_LWd = Fill_LWd(:,1); 
        elseif lat_inp_2 == 2
            % Grab all DOKIA data
            Fill_LWd = Fill_LWd(:,2); 
        end
        
    elseif lat_inp == 3 % Fill with LDAS data
        % Find similar LDAS LW  and Data_Fill times
        [C,IA,IB] = intersect(datenum_round_off(TIME_out,'minute'),datenum_round_off(TIME_30min(:,7),'minute')); 
        error('check line below is correct!')
        LW_LDAS = data_NLDAS(IB(1):IB(end),2); 
        
        % Fill with LDAS as second station
        [Fill_LWd_out] = data_fill_NW_edits(TIME_out_all(IA(1):IA(end),:), [Var_raw(IA(1):IA(end)) LW_LDAS], 1, 500, 100, 1, 1, 48, [site_E site_E]);
        % Grab only the first station (combined timeseries)
 
        
        % Plot check
        figure(1)
        hold on
        plot(TIME_out_all(:,7),Var_raw,'--b')
        plot(TIME_out_all(IA(1):IA(end),7),Fill_LWd_out(:,1),'b')
        plot(TIME_out_all(IA(1):IA(end),7),Fill_LWd_out(:,2),'k')
        legend('QCed Obs','Obs-LW-dw-FILLED','DOKIA-LW-dw')
        tlabel
        pause
        
        
        disp('Do want to use the filled data set, or ALL LDAS?')
        lat_inp_2 = user_inp(cellstr('Enter 1 - Filled. 2 - ALL LDAS'), 1, 2);
        if lat_inp_2 == 1
            % Grab only the first station (combined timeseries)
            Var_raw(IA(1):IA(end)) = Fill_LWd_out(:,1); 
            Fill_LWd = Var_raw;
        elseif lat_inp_2 == 2
            % Grab all LDAS data
            Var_raw(IA(1):IA(end)) = Fill_LWd_out(:,2); 
            Fill_LWd = Var_raw;
        end
        
        
    end
    
    disp('WARNING, need to add LW to data_fill_NW_edits')
    % Save Var Fill data
    eval(['save ' strcat(Filldir,Var_file_save) ' Fill_LWd']);
elseif u_inp == 1
   load(strcat(Filldir,Var_file_save)) 
end
Data_Fill(:,VAR_N) = Fill_LWd;
clc; clear alternate_VAR_N

%% END FORCING VARIABLES






%% VERIFICATION VARIABLES

%% MVF - 10 - Shortwave up
disp('Shortwave up')
VAR_N = 10;
Var_file_save = '\Fill_SWu.mat';

if exist(strcat(Filldir,Var_file_save)) == 2
    disp('Previous Filled Shortwave data found')
    u_inp = user_inp(cellstr('Enter 1 - Use this data. 2 - Create new Fill data. 3 - Continue with previous filled data'), 1, 2, 3);
else
    disp('No previous data found, starting filling procedure')
    u_inp = 2;
end

if u_inp == 2 || u_inp == 3
    if u_inp == 2
        % Get raw QC data
        Var_raw = Data_QC(:,VAR_N);
    elseif u_inp == 3 % Get previous Filled data
        load(strcat(Filldir,Var_file_save)) 
        Var_raw = Fill_SWu;
    end
    % Fill (Interactive)
    [Fill_SWu] = data_fill_NW_edits(TIME_out_all, Var_raw, 1, 1365, 0, 1, 1, 48, site_E);
    disp('WARNING, need to add SW to data_fill_NW_edits')
    % Save Var Fill data
    eval(['save ' strcat(Filldir,Var_file_save) ' Fill_SWu']);
elseif u_inp == 1
   load(strcat(Filldir,Var_file_save)) 
end
Data_Fill(:,VAR_N) = Fill_SWu;
clc; clear alternate_VAR_N


%% MVF - 12 - Longwave up
disp('Longwave up')
VAR_N = 12;
Var_file_save = '\Fill_LWu.mat';

if exist(strcat(Filldir,Var_file_save)) == 2
    disp('Previous Filled Longwave data found')
    u_inp = user_inp(cellstr('Enter 1 - Use this data. 2 - Create new Fill data. 3 - Continue with previous filled data'), 1, 2, 3);
else
    disp('No previous data found, starting filling procedure')
    u_inp = 2;
end

if u_inp == 2 || u_inp == 3
    if u_inp == 2
        % Get raw QC data
        Var_raw = Data_QC(:,VAR_N);
    elseif u_inp == 3 % Get previous Filled data
        load(strcat(Filldir,Var_file_save)) 
        Var_raw = Fill_LWu;
    end
    % Fill (Interactive)
    [Fill_LWu] = data_fill_NW_edits(TIME_out_all, Var_raw, 1, 500, 100, 1, 1, 48, site_E);
    disp('WARNING, need to add LW to data_fill_NW_edits')
    % Save Var Fill data
    eval(['save ' strcat(Filldir,Var_file_save) ' Fill_LWu']);
elseif u_inp == 1
   load(strcat(Filldir,Var_file_save)) 
end
Data_Fill(:,VAR_N) = Fill_LWu;
clc; clear alternate_VAR_N

%% MVF - 22 - Surface Temperature
disp('Surface Temperature')
VAR_N = 22;
Var_file_save = '\Fill_ST.mat';

if exist(strcat(Filldir,Var_file_save)) == 2
    disp('Previous Filled Surface Temperature data found')
    u_inp = user_inp(cellstr('Enter 1 - Use this data. 2 - Create new Fill data. 3 - Continue with previous filled data'), 1, 2, 3);
else
    disp('No previous data found, starting filling procedure')
    u_inp = 2;
end

if u_inp == 2 || u_inp == 3
    if u_inp == 2
        % Get raw QC data
        Var_raw = Data_QC(:,VAR_N);
    elseif u_inp == 3 % Get previous Filled data
        load(strcat(Filldir,Var_file_save)) 
        Var_raw = Fill_ST;
    end
    % Fill (Interactive)
    [Fill_ST] = data_fill_NW_edits(TIME_out_all, Var_raw, 1, 500, 100, 1, 1, 48, site_E);
    disp('WARNING, need to add LW to data_fill_NW_edits')
    % Save Var Fill data
    eval(['save ' strcat(Filldir,Var_file_save) ' Fill_ST']);
elseif u_inp == 1
   load(strcat(Filldir,Var_file_save)) 
end
Data_Fill(:,VAR_N) = Fill_ST;
clc; clear alternate_VAR_N

%% MVF - 20 - Snow Depth A m
disp('Snow Depth A m')
VAR_N = 20;
Var_file_save = '\Fill_SD.mat';

if exist(strcat(Filldir,Var_file_save)) == 2
    disp('Previous Filled Snow Depth m A data found')
    u_inp = user_inp(cellstr('Enter 1 - Use this data. 2 - Create new Fill data. 3 - Continue with previous filled data'), 1, 2, 3);
else
    disp('No previous data found, starting filling procedure')
    u_inp = 2;
end

if u_inp == 2 || u_inp == 3
    if u_inp == 2
        % Get raw QC data
        Var_raw = Data_QC(:,VAR_N);
    elseif u_inp == 3 % Get previous Filled data
        load(strcat(Filldir,Var_file_save)) 
        Var_raw = Fill_SD;
    end
    % Fill (Interactive)
    [Fill_SD] = data_fill_NW_edits(TIME_out_all, Var_raw, 1, 3, 0, 1, 1, 48, site_E);
    % Save Var Fill data
    eval(['save ' strcat(Filldir,Var_file_save) ' Fill_SD']);
elseif u_inp == 1
   load(strcat(Filldir,Var_file_save)) 
end
Data_Fill(:,VAR_N) = Fill_SD;
clc; clear alternate_VAR_N

%% MVF - 44 - Height of CSAT m
disp('Height of CSAT m')
VAR_N = 44;
Var_file_save = '\Fill_HCSAT.mat';

if exist(strcat(Filldir,Var_file_save)) == 2
    disp('Previous Filled Height of CSAT data found')
    u_inp = user_inp(cellstr('Enter 1 - Use this data. 2 - Create new Fill data. 3 - Continue with previous filled data'), 1, 2, 3);
else
    disp('No previous data found, starting filling procedure')
    u_inp = 2;
end

if u_inp == 2 || u_inp == 3
    if u_inp == 2
        % Get raw QC data
        Var_raw = Data_QC(:,VAR_N);
    elseif u_inp == 3 % Get previous Filled data
        load(strcat(Filldir,Var_file_save)) 
        Var_raw = Fill_HCSAT;
    end
    % Fill (Interactive)
    [Fill_HCSAT] = data_fill_NW_edits(TIME_out_all, Var_raw, 2, 10, 0, 1, 1, 48, site_E);
    disp('WARNING, need to add LW to data_fill_NW_edits')
    % Save Var Fill data
    eval(['save ' strcat(Filldir,Var_file_save) ' Fill_HCSAT']);
elseif u_inp == 1
   load(strcat(Filldir,Var_file_save)) 
end
Data_Fill(:,VAR_N) = Fill_HCSAT;
clc; clear alternate_VAR_N

%% MVF - 45 - CSAT Azimuth from true north
disp('CSAT Azimuth from true north')
VAR_N = 45;
Var_file_save = '\Fill_AzCSAT.mat';

if exist(strcat(Filldir,Var_file_save)) == 2
    disp('Previous CSAT Azimuth from true north data found')
    u_inp = user_inp(cellstr('Enter 1 - Use this data. 2 - Create new Fill data. 3 - Continue with previous filled data'), 1, 2, 3);
else
    disp('No previous data found, starting filling procedure')
    u_inp = 2;
end

if u_inp == 2 || u_inp == 3
    if u_inp == 2
        % Get raw QC data
        Var_raw = Data_QC(:,VAR_N);
    elseif u_inp == 3 % Get previous Filled data
        load(strcat(Filldir,Var_file_save)) 
        Var_raw = Fill_AzCSAT;
    end
    % Fill (Interactive)
    [Fill_AzCSAT] = data_fill_NW_edits(TIME_out_all, Var_raw, 5, 360, 0, 1, 1, 48, site_E);
    disp('WARNING, need to add LW to data_fill_NW_edits')
    % Save Var Fill data
    eval(['save ' strcat(Filldir,Var_file_save) ' Fill_AzCSAT']);
elseif u_inp == 1
   load(strcat(Filldir,Var_file_save)) 
end
Data_Fill(:,VAR_N) = Fill_AzCSAT;
clc; clear alternate_VAR_N


if Add_Spinup
    disp('............ Adding Spinup Period ................')
    
    %% Select Day to repeat for Spinup Period
    
    % Choose a day with no precip and to minimize melting of snow pack
    
    % Currently choose day manually
    
% % % %     figure(1)
% % % %     subplot(2,1,1)
% % % %     hold on
% % % %     plot(TIME_out,Data_Fill(:,1))
% % % %     r1 = refline(0,0);
% % % %     ylabel('Temperature (C)')
% % % %     tlabel
% % % %     
% % % %     subplot(2,1,2)
% % % %     hold on
% % % %     plot(TIME_out,Data_Fill(:,9))
% % % %     ylabel('SW down (W/m^2)')
% % % %     tlabel
    
    S_SP     = datenum(2012,12,8,1,0,0); % needs to be same hour/minute as first step in real data
    E_SP     = datenum(2012,12,9,1,0,0);
    
    IS_SP = find(TIME_out==S_SP); if isempty(IS_SP); error('not found'); end;
    IE_SP = find(TIME_out==E_SP); if isempty(IE_SP); error('not found'); end;
    
    NSUdays  = 30;
    tsperday = 24/dt_in;
        
    % Create Spinup timeseries
    TIME_SU = time_builder(TIME_out(1)-NSUdays,TIME_out(1)-((dt_in/2)/24),dt_in); % only subract a half of time step 
    SU_data = nan(size(TIME_SU,1),size(Data_Fill,2));
    I_mid = find(TIME_SU(:,7)==floor(TIME_SU(:,7))); 
    
    % Loop through each day in SU period, add SU day
    for cday = 1:length(I_mid)
        
        SU_data(  I_mid(cday) : I_mid(cday)+tsperday-1  ,:) = Data_Fill(IS_SP:IE_SP-1 , :); 
 
    end
    
    % Remove Preicpitation
    SU_data(:,2) = 0; 
    
    % Make winds Calm (0.5 m/s)
    SU_data(:,5) = 0.01; 
    
    % Make LW average over period
    SU_data(:,11) = mean(Data_Fill(IS_SP:IE_SP-1,11)); 
    
    %% Check EB consistency of Spin up period
    
    % Simple net EB estimation (we want it close to zero over day)
    % dQ = SWd*(1-albedo) + LWd - LWu
    dQ = SU_data(:,9) * (1-0.8) + SU_data(:,11) - (5.67*10^-8)*(272.0^4);
    
%     figure(2)
%     hold on
%     plot(TIME_SU(:,7),dQ,'b')
%     ylabel('Simple EB, (W m^-^2')
%     tlabel


    % Combine Spin up data and Real data, and times
    Data_Fill = [SU_data; Data_Fill];
    TIME_out  = [TIME_SU(:,7); TIME_out]; 
    TIME_out_all = time_builder(TIME_out); 
    Missing_I    = [logical(ones(size(SU_data))); Missing_I];
    
end

%% CHECK Forcing Variables for any NaN (Only Warn)
Nlength = size(Data_Fill,1); 
for Fvars = [1 2 4 5 16 16 9 11]
   PerNaN = sum(isnan(Data_Fill(:,Fvars))) ./ Nlength * 100;
   if PerNaN > 0
       sprintf('WARNING: There remain %3.1f percent NAN in variable %d',PerNaN,Fvars')
       Err_found = true;
   end
end

if exist('Err_found','var')
    disp('An error was found, see warning above')
end


%% Save Fill data
file_name_out = strcat(Filldir,'\',file_out);
eval(['save ' file_name_out ' Data_Fill TIME_out Missing_I MHeight_MVF site_E'])
fprintf('Finished Filling data for forcing %s \n',Forcing_N_out)
disp('Filled data Saved')
