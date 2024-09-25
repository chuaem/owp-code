%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mlr_metabolicDrivers.m
% This script creates multiple linear regression models to evaluate the influence
% of different environmental drivers on the estimated metabolic rates.
%
% AUTHOR:
% Emily Chua
%
% DATE:
% First created: 7/22/2024
% Last updated: 8/14/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;close all;clc

rootpath = 'G:\My Drive\Postdoc\Work\SMIIL\';

%==========================================================================
% Import data
%==========================================================================

%====Import final QC'd data & diel analysis results========================
site = 'gull';
cd([rootpath,'open-water-platform-data\',site,'\cleaned\final-qc'])
load([site,'-cleaned.mat']);
gull_params = finalQC;
cd([rootpath,'diel-method\matlab-results\final-qc\',site])
load('diel_res.mat')
gull_metab = diel_dtd;

site = 'north';
cd([rootpath,'open-water-platform-data\',site,'\cleaned\final-qc'])
load([site,'-cleaned.mat']);
north_params = finalQC;
cd([rootpath,'diel-method\matlab-results\final-qc\',site])
load('diel_res.mat')
north_metab = diel_dtd;

site = 'south';
cd([rootpath,'open-water-platform-data\',site,'\cleaned\final-qc'])
load([site,'-cleaned.mat']);
south_params = finalQC;
cd([rootpath,'diel-method\matlab-results\final-qc\',site])
load('diel_res.mat')
south_metab = diel_dtd;

clearvars finalQC diel_dtd diel_obs

% Find daily tidal range for each site
gull_range = groupsummary(gull_params,"datetime_utc","day","range","depth");
gull_range.day_datetime_utc = datetime(cellstr(gull_range.day_datetime_utc),'TimeZone','UTC');
gull_range = table2timetable(gull_range);
gull_range.GroupCount = [];
gull_range.Properties.VariableNames = {'tidal'};

north_range = groupsummary(north_params,"datetime_utc","day","range","depth");
north_range.day_datetime_utc = datetime(cellstr(north_range.day_datetime_utc),'TimeZone','UTC');
north_range = table2timetable(north_range);
north_range.GroupCount = [];
north_range.Properties.VariableNames = {'tidal'};

south_range = groupsummary(south_params,"datetime_utc","day","range","depth");
south_range.day_datetime_utc = datetime(cellstr(south_range.day_datetime_utc),'TimeZone','UTC');
south_range = table2timetable(south_range);
south_range.GroupCount = [];
south_range.Properties.VariableNames = {'tidal'};

%====Import windspeed and PAR data=========================================
cd([rootpath,'physical-data\final-dataset'])
load('windspeed.mat')

cd([rootpath,'physical-data\final-dataset'])
load('par.mat')

%====Calculate daily means for each site===================================
gull_dailyAvg = retime(gull_params,'daily','mean');
north_dailyAvg = retime(north_params,'daily','mean');
south_dailyAvg = retime(south_params,'daily','mean');

wspd_dailyAvg = retime(era5Dat,'daily','mean');
par_dailyAvg = retime(parDat,'daily','mean');

cd([rootpath,'figures\stats-analyses'])
%%
%====Create time tables for each site======================================
%----Gull------------------------------------------------------------------
dt2 = dateshift(gull_metab.daystart_dt,'start','day');
gull_metab.daystart_dt = dt2;
gull_TT = synchronize(gull_dailyAvg,gull_range,wspd_dailyAvg,par_dailyAvg,gull_metab);
gull_TT = removevars(gull_TT,{'deployment','daystart_dt','dayend_dt','daylength','R_hourly','P_hourly','R_daily','P_daily','datetime_local','Tair','light_lux'});
gull_TT(gull_TT.datetime_utc < datetime(2021,6,29,"TimeZone","UTC"),:) = [];
gull_TT(gull_TT.datetime_utc > datetime(2024,6,4,"TimeZone","UTC"),:) = [];

% Set anomalous GPP and ER values to NaN
anomER = find(gull_TT.ER > 0);
anomGPP = find(gull_TT.GPP < 0);
% gull_TT{[anomER;anomGPP],:} = NaN;
gull_TT([anomER;anomGPP],:) = [];

% Remove row if any explanatory or response variable is NaN
gull_TT = rmmissing(gull_TT,'DataVariables',{'tidal','salinity','temperature','chla','turbidity','wspd','GPP','ER','NEM'});

% Get month numbers
mo = month(gull_TT.datetime_utc);

% Define the seasons by month
isWinter = ismember(mo,[12 1 2]);
isSpring = ismember(mo,[3 4 5]);
isSummer = ismember(mo,[6 7 8]);
isFall = ismember(mo,[9 10 11]);

indWinter = find(isWinter);
indSpring = find(isSpring);
indSummer = find(isSummer);
indFall = find(isFall);

% Add a column for season number
gull_TT.season(indWinter) = 1;
gull_TT.season(indSpring) = 2;
gull_TT.season(indSummer) = 3;
gull_TT.season(indFall) = 4;
gull_TT.season = categorical(gull_TT.season);

% Add a column for site number
gull_TT.site = repelem(2,height(gull_TT))';
gull_TT.site = categorical(gull_TT.site);

%----North-----------------------------------------------------------------
dt2 = dateshift(north_metab.daystart_dt,'start','day');
north_metab.daystart_dt = dt2;
north_TT = synchronize(north_dailyAvg,north_range,wspd_dailyAvg,par_dailyAvg,north_metab);
north_TT = removevars(north_TT,{'deployment','daystart_dt','dayend_dt','daylength','R_hourly','P_hourly','R_daily','P_daily','datetime_local','Tair','light_lux'});
north_TT(north_TT.datetime_utc < datetime(2021,6,29,"TimeZone","UTC"),:) = [];
north_TT(north_TT.datetime_utc > datetime(2024,6,4,"TimeZone","UTC"),:) = [];

% Set anomalous GPP and ER values to NaN
anomER = find(north_TT.ER > 0);
anomGPP = find(north_TT.GPP < 0);
% north_TT{[anomER;anomGPP],:} = NaN;
north_TT([anomER;anomGPP],:) = [];

% Remove row if any explanatory or response variable is NaN
north_TT = rmmissing(north_TT,'DataVariables',{'tidal','salinity','temperature','chla','turbidity','wspd','GPP','ER','NEM'});

% Get month numbers
mo = month(north_TT.datetime_utc);

% Define the seasons by month
isWinter = ismember(mo,[12 1 2]);
isSpring = ismember(mo,[3 4 5]);
isSummer = ismember(mo,[6 7 8]);
isFall = ismember(mo,[9 10 11]);

indWinter = find(isWinter);
indSpring = find(isSpring);
indSummer = find(isSummer);
indFall = find(isFall);

% Add a column for season number
north_TT.season(indWinter) = 1;
north_TT.season(indSpring) = 2;
north_TT.season(indSummer) = 3;
north_TT.season(indFall) = 4;
north_TT.season = categorical(north_TT.season);

% Add a column for site number
north_TT.site = repelem(1,height(north_TT))';
north_TT.site = categorical(north_TT.site);

%----South-----------------------------------------------------------------
dt2 = dateshift(south_metab.daystart_dt,'start','day');
south_metab.daystart_dt = dt2;
south_TT = synchronize(south_dailyAvg,south_range,wspd_dailyAvg,par_dailyAvg,south_metab);
south_TT = removevars(south_TT,{'deployment','daystart_dt','dayend_dt','daylength','R_hourly','P_hourly','R_daily','P_daily','datetime_local','Tair','light_lux'});
south_TT(south_TT.datetime_utc < datetime(2021,6,29,"TimeZone","UTC"),:) = [];
south_TT(south_TT.datetime_utc > datetime(2024,6,4,"TimeZone","UTC"),:) = [];

% Set anomalous GPP and ER values to NaN
anomER = find(south_TT.ER > 0);
anomGPP = find(south_TT.GPP < 0);
% south_TT{[anomER;anomGPP],:} = NaN;
south_TT([anomER;anomGPP],:) = [];

% Remove row if any explanatory or response variable is NaN
south_TT = rmmissing(south_TT,'DataVariables',{'tidal','salinity','temperature','chla','turbidity','wspd','GPP','ER','NEM'});

% Get month numbers
mo = month(south_TT.datetime_utc);

% Define the seasons by month
isWinter = ismember(mo,[12 1 2]);
isSpring = ismember(mo,[3 4 5]);
isSummer = ismember(mo,[6 7 8]);
isFall = ismember(mo,[9 10 11]);

indWinter = find(isWinter);
indSpring = find(isSpring);
indSummer = find(isSummer);
indFall = find(isFall);

% Add a column for season number
south_TT.season(indWinter) = 1;
south_TT.season(indSpring) = 2;
south_TT.season(indSummer) = 3;
south_TT.season(indFall) = 4;
south_TT.season = categorical(south_TT.season);

% Add a column for site number
south_TT.site = repelem(3,height(south_TT))';
south_TT.site = categorical(south_TT.site);

%====Create a table to store predictor variables for all sites=============
gull_dat.all = timetable2table(gull_TT);
gull_dat.all.datetime_utc = [];

north_dat.all = timetable2table(north_TT);
north_dat.all.datetime_utc = [];

south_dat.all = timetable2table(south_TT);
south_dat.all.datetime_utc = [];

allSites_dat.all = [gull_dat.all;north_dat.all;south_dat.all];

% % HERE
% idx = find(allSites_dat.all.NEM < -60);
% allSites_dat.all(idx,:) = [];

%==========================================================================
% Construct multiple linear regression models for ER, GPP, and NEM
%==========================================================================

%====MLRs across all sites: (1) for each season and (2) across all seasons=

%----Original (not scaled) data--------------------------------------------
% Break data into seasons
isWinter = ismember(allSites_dat.all.season,'1');
isSpring = ismember(allSites_dat.all.season,'2');
isSummer = ismember(allSites_dat.all.season,'3');
isFall = ismember(allSites_dat.all.season,'4');

indWinter = find(isWinter);
indSpring = find(isSpring);
indSummer = find(isSummer);
indFall = find(isFall);

allSites_dat.winter = allSites_dat.all(indWinter,:);
allSites_dat.spring = allSites_dat.all(indSpring,:);
allSites_dat.summer = allSites_dat.all(indSummer,:);
allSites_dat.fall = allSites_dat.all(indFall,:);

% Compute models with original data
data = allSites_dat;

% (1) By season
% GPP
allSites_mdl.winter.GPP = fitlm(data.winter,'GPP ~ tidal + wspd + temperature + salinity + chla + turbidity + site');
allSites_mdl.spring.GPP = fitlm(data.spring,'GPP ~ tidal + wspd + temperature + salinity + chla + turbidity + site');
allSites_mdl.summer.GPP = fitlm(data.summer,'GPP ~ tidal + wspd + temperature + salinity + chla + turbidity + site');
allSites_mdl.fall.GPP = fitlm(data.fall,'GPP ~ tidal + wspd + temperature + salinity + chla + turbidity + site');
% ER
allSites_mdl.winter.ER = fitlm(data.winter,'ER ~ tidal + wspd + temperature + salinity + chla + turbidity + site');
allSites_mdl.spring.ER = fitlm(data.spring,'ER ~ tidal + wspd + temperature + salinity + chla + turbidity + site');
allSites_mdl.summer.ER = fitlm(data.summer,'ER ~ tidal + wspd + temperature + salinity + chla + turbidity + site');
allSites_mdl.fall.ER = fitlm(data.fall,'ER ~ tidal + wspd + temperature + salinity + chla + turbidity + site');
% NEM
allSites_mdl.winter.NEM = fitlm(data.winter,'NEM ~ tidal + wspd + temperature + salinity + chla + turbidity + site');
allSites_mdl.spring.NEM = fitlm(data.spring,'NEM ~ tidal + wspd + temperature + salinity + chla + turbidity + site');
allSites_mdl.summer.NEM = fitlm(data.summer,'NEM ~ tidal + wspd + temperature + salinity + chla + turbidity + site');
allSites_mdl.fall.NEM = fitlm(data.fall,'NEM ~ tidal + wspd + temperature + salinity + chla + turbidity + site');

% (2) Across all seasons
% GPP
allSites_mdl.all.GPP = fitlm(data.all,'GPP ~ tidal + wspd + temperature + salinity + chla + turbidity + season*site'); % Include interactions
% ER
allSites_mdl.all.ER = fitlm(data.all,'ER ~ tidal + wspd + temperature + salinity + chla + turbidity + season*site'); % Include interactions
% NEM
allSites_mdl.all.NEM = fitlm(data.all,'NEM ~ tidal + wspd + temperature + salinity + chla + turbidity + season*site'); % Include interactions

%----Scaled data-----------------------------------------------------------
% Scale response and explanatory variables to center 0 and std 1 (Lowe et al., 2019)
VN = data.all.Properties.VariableNames;
colNr = find(strcmp(VN,'NEM'));   % Find column number for NEM
allSites_dat_norm.all = normalize(allSites_dat.all(:,1:colNr)); % Normalize all parameters except "season" and "site"
allSites_dat_norm.all.("season") = data.all.season; % Add "season" column back in
allSites_dat_norm.all.("site") = data.all.site; % Add "site" column back in

allSites_dat_norm.winter = allSites_dat_norm.all(indWinter,:);
allSites_dat_norm.spring = allSites_dat_norm.all(indSpring,:);
allSites_dat_norm.summer = allSites_dat_norm.all(indSummer,:);
allSites_dat_norm.fall = allSites_dat_norm.all(indFall,:);

% Compute models with normalized data
data_norm = allSites_dat_norm;

% (1) By season
% GPP
allSites_mdl_norm.winter.GPP = fitlm(data_norm.winter,'GPP ~ tidal + wspd + temperature + salinity + chla + turbidity + site');
allSites_mdl_norm.spring.GPP = fitlm(data_norm.spring,'GPP ~ tidal + wspd + temperature + salinity + chla + turbidity + site');
allSites_mdl_norm.summer.GPP = fitlm(data_norm.summer,'GPP ~ tidal + wspd + temperature + salinity + chla + turbidity + site');
allSites_mdl_norm.fall.GPP = fitlm(data_norm.fall,'GPP ~ tidal + wspd + temperature + salinity + chla + turbidity + site');
% ER
allSites_mdl_norm.winter.ER = fitlm(data_norm.winter,'ER ~ tidal + wspd + temperature + salinity + chla + turbidity + site');
allSites_mdl_norm.spring.ER = fitlm(data_norm.spring,'ER ~ tidal + wspd + temperature + salinity + chla + turbidity + site');
allSites_mdl_norm.summer.ER = fitlm(data_norm.summer,'ER ~ tidal + wspd + temperature + salinity + chla + turbidity + site');
allSites_mdl_norm.fall.ER = fitlm(data_norm.fall,'ER ~ tidal + wspd + temperature + salinity + chla + turbidity + site');
% NEM
allSites_mdl_norm.winter.NEM = fitlm(data_norm.winter,'NEM ~ tidal + wspd + temperature + salinity + chla + turbidity + site');
allSites_mdl_norm.spring.NEM = fitlm(data_norm.spring,'NEM ~ tidal + wspd + temperature + salinity + chla + turbidity + site');
allSites_mdl_norm.summer.NEM = fitlm(data_norm.summer,'NEM ~ tidal + wspd + temperature + salinity + chla + turbidity + site');
allSites_mdl_norm.fall.NEM = fitlm(data_norm.fall,'NEM ~ tidal + wspd + temperature + salinity + chla + turbidity + site');

% (2) Across all seasons
% GPP
allSites_mdl_norm.all.GPP = fitlm(data_norm.all,'GPP ~ tidal + wspd + temperature + salinity + chla + turbidity + season*site'); % Include interactions
% ER
allSites_mdl_norm.all.ER = fitlm(data_norm.all,'ER ~ tidal + wspd + temperature + salinity + chla + turbidity + season*site'); % Include interactions
% NEM
allSites_mdl_norm.all.NEM = fitlm(data_norm.all,'NEM ~ tidal + wspd + temperature + salinity + chla + turbidity + season*site'); % Include interactions

% % MIXED EFFECTS MODEL
% allSites_mdl_norm_MIXED.all.GPP = fitlme(data_norm.all,'GPP ~ tidal + wspd + temperature + salinity + chla + turbidity + site*season'); % Gives same thing as fitlm
% 
% allSites_mdl_norm_MIXED2.all.GPP = fitlme(data_norm.all,'GPP ~ tidal + wspd + temperature + salinity + chla + turbidity + season*site + (site|season)');

%%
%==========================================================================
% Plot results
%==========================================================================
data = allSites_dat;
mdl = allSites_mdl;

respNames = {'GPP','ER','NEM'};
predNames = {'salinity','temperature','chla','turbidity','tidal','wspd'};
seasonNames = {'all','winter','spring','summer','fall'};
xlblNames = {'Salinity','Temperature (^oC)','Chl a (RFU)','Turbidity (NTU)','Tidal (m)','Wind speed (m/s)'};

% Find means of each variable for plotting individual regression lines
for k = 1:length(predNames)
    meanVal.all.(predNames{k}) = mean(data.all.(predNames{k}),'omitmissing');
    meanVal.winter.(predNames{k}) = mean(data.winter.(predNames{k}),'omitmissing');
    meanVal.spring.(predNames{k}) = mean(data.spring.(predNames{k}),'omitmissing');
    meanVal.summer.(predNames{k}) = mean(data.summer.(predNames{k}),'omitmissing');
    meanVal.fall.(predNames{k}) = mean(data.fall.(predNames{k}),'omitmissing');
end
%%
%====GPP, ER, NEM vs. Explanatory Variables================================
for i = 1:length(respNames)   % Loop through plots for GPP, ER, and NEM
    fig(i) = figure(i);clf
    % fig(i).Position = [363.4 905.0 676.8 957.6];
    % fig(i).Position = [363.4 875.4 655.2 957.6];
    fig(i).Position = [-1337.4 -189.4 664 957.6];
    t = tiledlayout(3,2,'Padding','none','TileSpacing','loose');

    % Coefficient estimates for each "season"
    for m = 1:length(seasonNames)
        % Model with all data
        b0.(seasonNames{m}) = mdl.(seasonNames{m}).(respNames{i}).Coefficients.Estimate(1);
        b1.(seasonNames{m}) = mdl.(seasonNames{m}).(respNames{i}).Coefficients.Estimate(2);
        b2.(seasonNames{m}) = mdl.(seasonNames{m}).(respNames{i}).Coefficients.Estimate(3);
        b3.(seasonNames{m}) = mdl.(seasonNames{m}).(respNames{i}).Coefficients.Estimate(4);
        b4.(seasonNames{m}) = mdl.(seasonNames{m}).(respNames{i}).Coefficients.Estimate(5);
        b5.(seasonNames{m}) = mdl.(seasonNames{m}).(respNames{i}).Coefficients.Estimate(6);
        b6.(seasonNames{m}) = mdl.(seasonNames{m}).(respNames{i}).Coefficients.Estimate(7);
    end

    for j = 1:length(predNames) % Loop through plotting each predictor
        t1 = tiledlayout(t,2,2,'Tilespacing','none');
        t1.Layout.Tile = j;

        x = data.all.(predNames{j});

        % Winter
        ax1 = nexttile(t1);
        plot(data.all.(predNames{j}),data.all.(respNames{i}),'.','color',rgb('lightgrey'))   % Plot All data in grey
        hold on
        % Plot regression line for All data in grey
        switch predNames{j}
            case 'salinity'
                plot(x, b0.all + mdl.all.(respNames{i}).Coefficients.Estimate(j+1)*x...
                    + b2.all*meanVal.all.(predNames{2}) + b3.all*meanVal.all.(predNames{3}) + b4.all*meanVal.all.(predNames{4})...
                    + b5.all*meanVal.all.(predNames{5}) + b6.all*meanVal.all.(predNames{6}),'color',rgb('grey'))
            case 'temperature'
                plot(x, b0.all + b1.all*meanVal.all.(predNames{1}) + mdl.all.(respNames{i}).Coefficients.Estimate(j+1)*x...
                    + b3.all*meanVal.all.(predNames{3}) + b4.all*meanVal.all.(predNames{4})...
                    + b5.all*meanVal.all.(predNames{5}) + b6.all*meanVal.all.(predNames{6}),'color',rgb('grey'))
            case 'chla'
                plot(x, b0.all + b1.all*meanVal.all.(predNames{1}) + b2.all*meanVal.all.(predNames{2})...
                    + mdl.all.(respNames{i}).Coefficients.Estimate(j+1)*x + b4.all*meanVal.all.(predNames{4})...
                    + b5.all*meanVal.all.(predNames{5}) + b6.all*meanVal.all.(predNames{6}),'color',rgb('grey'))
            case 'turbidity'
                plot(x, b0.all + b1.all*meanVal.all.(predNames{1}) + b2.all*meanVal.all.(predNames{2})...
                    + b3.all*meanVal.all.(predNames{3}) + mdl.all.(respNames{i}).Coefficients.Estimate(j+1)*x...
                    + b5.all*meanVal.all.(predNames{5}) + b6.all*meanVal.all.(predNames{6}),'color',rgb('grey'))
            case 'tidal'
                plot(x, b0.all + b1.all*meanVal.all.(predNames{1}) + b2.all*meanVal.all.(predNames{2})...
                    + b3.all*meanVal.all.(predNames{3}) + b4.all*meanVal.all.(predNames{4})...
                    + mdl.all.(respNames{i}).Coefficients.Estimate(j+1)*x + b6.all*meanVal.all.(predNames{6}),'color',rgb('grey'))
            case 'wspd'
                plot(x, b0.all + b1.all*meanVal.all.(predNames{1}) + b2.all*meanVal.all.(predNames{2})...
                    + b3.all*meanVal.all.(predNames{3}) + b4.all*meanVal.all.(predNames{4})...
                    + b5.all*meanVal.all.(predNames{5}) + mdl.all.(respNames{i}).Coefficients.Estimate(j+1)*x,'color',rgb('grey'))
        end
        % Plot Winter data in black
        plot(data.winter.(predNames{j}),data.winter.(respNames{i}),'.k')
        % Plot regression line for Winter data in black
        switch predNames{j}
            case 'salinity'
                plot(x, b0.winter + mdl.winter.(respNames{i}).Coefficients.Estimate(j+1)*x...
                    + b2.winter*meanVal.winter.(predNames{2}) + b3.winter*meanVal.winter.(predNames{3}) + b4.winter*meanVal.winter.(predNames{4})...
                    + b5.winter*meanVal.winter.(predNames{5}),'k')
            case 'temperature'
                plot(x, b0.winter + b1.winter*meanVal.winter.(predNames{1}) + mdl.winter.(respNames{i}).Coefficients.Estimate(j+1)*x...
                    + b3.winter*meanVal.winter.(predNames{3}) + b4.winter*meanVal.winter.(predNames{4})...
                    + b5.winter*meanVal.winter.(predNames{5}),'k')
            case 'chla'
                plot(x, b0.winter + b1.winter*meanVal.winter.(predNames{1}) + b2.winter*meanVal.winter.(predNames{2})...
                    + mdl.winter.(respNames{i}).Coefficients.Estimate(j+1)*x + b4.winter*meanVal.winter.(predNames{4})...
                    + b5.winter*meanVal.winter.(predNames{5}) + b6.winter*meanVal.winter.(predNames{6}),'k')
            case 'turbidity'
                plot(x, b0.winter + b1.winter*meanVal.winter.(predNames{1}) + b2.winter*meanVal.winter.(predNames{2})...
                    + b3.winter*meanVal.winter.(predNames{3}) + mdl.winter.(respNames{i}).Coefficients.Estimate(j+1)*x...
                    + b5.winter*meanVal.winter.(predNames{5}) + b6.winter*meanVal.winter.(predNames{6}),'k')
            case 'tidal'
                plot(x, b0.winter + b1.winter*meanVal.winter.(predNames{1}) + b2.winter*meanVal.winter.(predNames{2})...
                    + b3.winter*meanVal.winter.(predNames{3}) + b4.winter*meanVal.winter.(predNames{4})...
                    + mdl.winter.(respNames{i}).Coefficients.Estimate(j+1)*x + b6.winter*meanVal.winter.(predNames{6}),'k')
            case 'wspd'
                plot(x, b0.winter + b1.winter*meanVal.winter.(predNames{1}) + b2.winter*meanVal.winter.(predNames{2})...
                    + b3.winter*meanVal.winter.(predNames{3}) + b4.winter*meanVal.winter.(predNames{4})...
                    + b5.winter*meanVal.winter.(predNames{5}) + mdl.winter.(respNames{i}).Coefficients.Estimate(j+1)*x,'k')
        end
        pbaspect(ax1,[1 1 1]); % Make relative lengths of axes equal
        % make a text object for the title
        xl = get(gca(),'Xlim');
        yl = get(gca(),'Ylim');
        text(xl(1),yl(2),'Winter','VerticalAlignment','top','HorizontalAlignment','left','FontSize',12)
        % text(xl(1)*1.05,yl(2)*0.95,'Winter','VerticalAlignment','top','HorizontalAlignment','left','FontSize',12)

        % Spring
        ax2 = nexttile(t1);
        plot(data.all.(predNames{j}),data.all.(respNames{i}),'.','color',rgb('lightgrey'))
        hold on
        % Plot regression line for All data in grey
        switch predNames{j}
            case 'salinity'
                plot(x, b0.all + mdl.all.(respNames{i}).Coefficients.Estimate(j+1)*x...
                    + b2.all*meanVal.all.(predNames{2}) + b3.all*meanVal.all.(predNames{3}) + b4.all*meanVal.all.(predNames{4})...
                    + b5.all*meanVal.all.(predNames{5}) + b6.all*meanVal.all.(predNames{6}),'color',rgb('grey'))
            case 'temperature'
                plot(x, b0.all + b1.all*meanVal.all.(predNames{1}) + mdl.all.(respNames{i}).Coefficients.Estimate(j+1)*x...
                    + b3.all*meanVal.all.(predNames{3}) + b4.all*meanVal.all.(predNames{4})...
                    + b5.all*meanVal.all.(predNames{5}) + b6.all*meanVal.all.(predNames{6}),'color',rgb('grey'))
            case 'chla'
                plot(x, b0.all + b1.all*meanVal.all.(predNames{1}) + b2.all*meanVal.all.(predNames{2})...
                    + mdl.all.(respNames{i}).Coefficients.Estimate(j+1)*x + b4.all*meanVal.all.(predNames{4})...
                    + b5.all*meanVal.all.(predNames{5}) + b6.all*meanVal.all.(predNames{6}),'color',rgb('grey'))
            case 'turbidity'
                plot(x, b0.all + b1.all*meanVal.all.(predNames{1}) + b2.all*meanVal.all.(predNames{2})...
                    + b3.all*meanVal.all.(predNames{3}) + mdl.all.(respNames{i}).Coefficients.Estimate(j+1)*x...
                    + b5.all*meanVal.all.(predNames{5}) + b6.all*meanVal.all.(predNames{6}),'color',rgb('grey'))
            case 'tidal'
                plot(x, b0.all + b1.all*meanVal.all.(predNames{1}) + b2.all*meanVal.all.(predNames{2})...
                    + b3.all*meanVal.all.(predNames{3}) + b4.all*meanVal.all.(predNames{4})...
                    + mdl.all.(respNames{i}).Coefficients.Estimate(j+1)*x + b6.all*meanVal.all.(predNames{6}),'color',rgb('grey'))
            case 'wspd'
                plot(x, b0.all + b1.all*meanVal.all.(predNames{1}) + b2.all*meanVal.all.(predNames{2})...
                    + b3.all*meanVal.all.(predNames{3}) + b4.all*meanVal.all.(predNames{4})...
                    + b5.all*meanVal.all.(predNames{5}) + mdl.all.(respNames{i}).Coefficients.Estimate(j+1)*x,'color',rgb('grey'))
        end
        % Plot Spring data in black
        plot(data.spring.(predNames{j}),data.spring.(respNames{i}),'.k')
        % Plot regression line for Spring data in black
        switch predNames{j}
            case 'salinity'
                plot(x, b0.spring + mdl.spring.(respNames{i}).Coefficients.Estimate(j+1)*x...
                    + b2.spring*meanVal.spring.(predNames{2}) + b3.spring*meanVal.spring.(predNames{3}) + b4.spring*meanVal.spring.(predNames{4})...
                    + b5.spring*meanVal.spring.(predNames{5}),'k')
            case 'temperature'
                plot(x, b0.spring + b1.spring*meanVal.spring.(predNames{1}) + mdl.spring.(respNames{i}).Coefficients.Estimate(j+1)*x...
                    + b3.spring*meanVal.spring.(predNames{3}) + b4.spring*meanVal.spring.(predNames{4})...
                    + b5.spring*meanVal.spring.(predNames{5}),'k')
            case 'chla'
                plot(x, b0.spring + b1.spring*meanVal.spring.(predNames{1}) + b2.spring*meanVal.spring.(predNames{2})...
                    + mdl.spring.(respNames{i}).Coefficients.Estimate(j+1)*x + b4.spring*meanVal.spring.(predNames{4})...
                    + b5.spring*meanVal.spring.(predNames{5}) + b6.spring*meanVal.spring.(predNames{6}),'k')
            case 'turbidity'
                plot(x, b0.spring + b1.spring*meanVal.spring.(predNames{1}) + b2.spring*meanVal.spring.(predNames{2})...
                    + b3.spring*meanVal.spring.(predNames{3}) + mdl.spring.(respNames{i}).Coefficients.Estimate(j+1)*x...
                    + b5.spring*meanVal.spring.(predNames{5}) + b6.spring*meanVal.spring.(predNames{6}),'k')
            case 'tidal'
                plot(x, b0.spring + b1.spring*meanVal.spring.(predNames{1}) + b2.spring*meanVal.spring.(predNames{2})...
                    + b3.spring*meanVal.spring.(predNames{3}) + b4.spring*meanVal.spring.(predNames{4})...
                    + mdl.spring.(respNames{i}).Coefficients.Estimate(j+1)*x + b6.spring*meanVal.spring.(predNames{6}),'k')
            case 'wspd'
                plot(x, b0.spring + b1.spring*meanVal.spring.(predNames{1}) + b2.spring*meanVal.spring.(predNames{2})...
                    + b3.spring*meanVal.spring.(predNames{3}) + b4.spring*meanVal.spring.(predNames{4})...
                    + b5.spring*meanVal.spring.(predNames{5}) + mdl.spring.(respNames{i}).Coefficients.Estimate(j+1)*x,'k')
        end
        pbaspect(ax2,[1 1 1]); % Make relative lengths of axes equal
        % make a text object for the title
        xl = get(gca(),'Xlim');
        yl = get(gca(),'Ylim');
        text(xl(1),yl(2),'Spring','VerticalAlignment','top','HorizontalAlignment','left','FontSize',12)
        % text(xl(1)*1.05,yl(2)*0.95,'Spring','VerticalAlignment','top','HorizontalAlignment','left','FontSize',12)

        % Summer
        ax3 = nexttile(t1);
        plot(data.all.(predNames{j}),data.all.(respNames{i}),'.','color',rgb('lightgrey'))
        hold on
        % Plot regression line for All data in grey
        switch predNames{j}
            case 'salinity'
                plot(x, b0.all + mdl.all.(respNames{i}).Coefficients.Estimate(j+1)*x...
                    + b2.all*meanVal.all.(predNames{2}) + b3.all*meanVal.all.(predNames{3}) + b4.all*meanVal.all.(predNames{4})...
                    + b5.all*meanVal.all.(predNames{5}) + b6.all*meanVal.all.(predNames{6}),'color',rgb('grey'))
            case 'temperature'
                plot(x, b0.all + b1.all*meanVal.all.(predNames{1}) + mdl.all.(respNames{i}).Coefficients.Estimate(j+1)*x...
                    + b3.all*meanVal.all.(predNames{3}) + b4.all*meanVal.all.(predNames{4})...
                    + b5.all*meanVal.all.(predNames{5}) + b6.all*meanVal.all.(predNames{6}),'color',rgb('grey'))
            case 'chla'
                plot(x, b0.all + b1.all*meanVal.all.(predNames{1}) + b2.all*meanVal.all.(predNames{2})...
                    + mdl.all.(respNames{i}).Coefficients.Estimate(j+1)*x + b4.all*meanVal.all.(predNames{4})...
                    + b5.all*meanVal.all.(predNames{5}) + b6.all*meanVal.all.(predNames{6}),'color',rgb('grey'))
            case 'turbidity'
                plot(x, b0.all + b1.all*meanVal.all.(predNames{1}) + b2.all*meanVal.all.(predNames{2})...
                    + b3.all*meanVal.all.(predNames{3}) + mdl.all.(respNames{i}).Coefficients.Estimate(j+1)*x...
                    + b5.all*meanVal.all.(predNames{5}) + b6.all*meanVal.all.(predNames{6}),'color',rgb('grey'))
            case 'tidal'
                plot(x, b0.all + b1.all*meanVal.all.(predNames{1}) + b2.all*meanVal.all.(predNames{2})...
                    + b3.all*meanVal.all.(predNames{3}) + b4.all*meanVal.all.(predNames{4})...
                    + mdl.all.(respNames{i}).Coefficients.Estimate(j+1)*x + b6.all*meanVal.all.(predNames{6}),'color',rgb('grey'))
            case 'wspd'
                plot(x, b0.all + b1.all*meanVal.all.(predNames{1}) + b2.all*meanVal.all.(predNames{2})...
                    + b3.all*meanVal.all.(predNames{3}) + b4.all*meanVal.all.(predNames{4})...
                    + b5.all*meanVal.all.(predNames{5}) + mdl.all.(respNames{i}).Coefficients.Estimate(j+1)*x,'color',rgb('grey'))
        end
        % Plot Summer data in black
        plot(data.summer.(predNames{j}),data.summer.(respNames{i}),'.k')
        % Plot regression line for Summer data in black
        switch predNames{j}
            case 'salinity'
                plot(x, b0.summer + mdl.summer.(respNames{i}).Coefficients.Estimate(j+1)*x...
                    + b2.summer*meanVal.summer.(predNames{2}) + b3.summer*meanVal.summer.(predNames{3}) + b4.summer*meanVal.summer.(predNames{4})...
                    + b5.summer*meanVal.summer.(predNames{5}),'k')
            case 'temperature'
                plot(x, b0.summer + b1.summer*meanVal.summer.(predNames{1}) + mdl.summer.(respNames{i}).Coefficients.Estimate(j+1)*x...
                    + b3.summer*meanVal.summer.(predNames{3}) + b4.summer*meanVal.summer.(predNames{4})...
                    + b5.summer*meanVal.summer.(predNames{5}),'k')
            case 'chla'
                plot(x, b0.summer + b1.summer*meanVal.summer.(predNames{1}) + b2.summer*meanVal.summer.(predNames{2})...
                    + mdl.summer.(respNames{i}).Coefficients.Estimate(j+1)*x + b4.summer*meanVal.summer.(predNames{4})...
                    + b5.summer*meanVal.summer.(predNames{5}) + b6.summer*meanVal.summer.(predNames{6}),'k')
            case 'turbidity'
                plot(x, b0.summer + b1.summer*meanVal.summer.(predNames{1}) + b2.summer*meanVal.summer.(predNames{2})...
                    + b3.summer*meanVal.summer.(predNames{3}) + mdl.summer.(respNames{i}).Coefficients.Estimate(j+1)*x...
                    + b5.summer*meanVal.summer.(predNames{5}) + b6.summer*meanVal.summer.(predNames{6}),'k')
            case 'tidal'
                plot(x, b0.summer + b1.summer*meanVal.summer.(predNames{1}) + b2.summer*meanVal.summer.(predNames{2})...
                    + b3.summer*meanVal.summer.(predNames{3}) + b4.summer*meanVal.summer.(predNames{4})...
                    + mdl.summer.(respNames{i}).Coefficients.Estimate(j+1)*x + b6.summer*meanVal.summer.(predNames{6}),'k')
            case 'wspd'
                plot(x, b0.summer + b1.summer*meanVal.summer.(predNames{1}) + b2.summer*meanVal.summer.(predNames{2})...
                    + b3.summer*meanVal.summer.(predNames{3}) + b4.summer*meanVal.summer.(predNames{4})...
                    + b5.summer*meanVal.summer.(predNames{5}) + mdl.summer.(respNames{i}).Coefficients.Estimate(j+1)*x,'k')
        end
        pbaspect(ax3,[1 1 1]); % Make relative lengths of axes equal
        % make a text object for the title
        xl = get(gca(),'Xlim');
        yl = get(gca(),'Ylim');
        text(xl(1),yl(2),'Summer','VerticalAlignment','top','HorizontalAlignment','left','FontSize',12)
        % text(xl(1)+.5,yl(2)-5,'Summer','VerticalAlignment','top','HorizontalAlignment','left','FontSize',12)

        % Fall
        ax4 = nexttile(t1);
        plot(data.all.(predNames{j}),data.all.(respNames{i}),'.','color',rgb('lightgrey'))
        hold on
        % Plot regression line for All data in grey
        switch predNames{j}
            case 'salinity'
                plot(x, b0.all + mdl.all.(respNames{i}).Coefficients.Estimate(j+1)*x...
                    + b2.all*meanVal.all.(predNames{2}) + b3.all*meanVal.all.(predNames{3}) + b4.all*meanVal.all.(predNames{4})...
                    + b5.all*meanVal.all.(predNames{5}) + b6.all*meanVal.all.(predNames{6}),'color',rgb('grey'))
            case 'temperature'
                plot(x, b0.all + b1.all*meanVal.all.(predNames{1}) + mdl.all.(respNames{i}).Coefficients.Estimate(j+1)*x...
                    + b3.all*meanVal.all.(predNames{3}) + b4.all*meanVal.all.(predNames{4})...
                    + b5.all*meanVal.all.(predNames{5}) + b6.all*meanVal.all.(predNames{6}),'color',rgb('grey'))
            case 'chla'
                plot(x, b0.all + b1.all*meanVal.all.(predNames{1}) + b2.all*meanVal.all.(predNames{2})...
                    + mdl.all.(respNames{i}).Coefficients.Estimate(j+1)*x + b4.all*meanVal.all.(predNames{4})...
                    + b5.all*meanVal.all.(predNames{5}) + b6.all*meanVal.all.(predNames{6}),'color',rgb('grey'))
            case 'turbidity'
                plot(x, b0.all + b1.all*meanVal.all.(predNames{1}) + b2.all*meanVal.all.(predNames{2})...
                    + b3.all*meanVal.all.(predNames{3}) + mdl.all.(respNames{i}).Coefficients.Estimate(j+1)*x...
                    + b5.all*meanVal.all.(predNames{5}) + b6.all*meanVal.all.(predNames{6}),'color',rgb('grey'))
            case 'tidal'
                plot(x, b0.all + b1.all*meanVal.all.(predNames{1}) + b2.all*meanVal.all.(predNames{2})...
                    + b3.all*meanVal.all.(predNames{3}) + b4.all*meanVal.all.(predNames{4})...
                    + mdl.all.(respNames{i}).Coefficients.Estimate(j+1)*x + b6.all*meanVal.all.(predNames{6}),'color',rgb('grey'))
            case 'wspd'
                plot(x, b0.all + b1.all*meanVal.all.(predNames{1}) + b2.all*meanVal.all.(predNames{2})...
                    + b3.all*meanVal.all.(predNames{3}) + b4.all*meanVal.all.(predNames{4})...
                    + b5.all*meanVal.all.(predNames{5}) + mdl.all.(respNames{i}).Coefficients.Estimate(j+1)*x,'color',rgb('grey'))
        end
        % Plot Fall data in black
        plot(data.fall.(predNames{j}),data.fall.(respNames{i}),'.k')
        % Plot regression line for Fall data in black
        switch predNames{j}
            case 'salinity'
                plot(x, b0.fall + mdl.fall.(respNames{i}).Coefficients.Estimate(j+1)*x...
                    + b2.fall*meanVal.fall.(predNames{2}) + b3.fall*meanVal.fall.(predNames{3}) + b4.fall*meanVal.fall.(predNames{4})...
                    + b5.fall*meanVal.fall.(predNames{5}),'k')
            case 'temperature'
                plot(x, b0.fall + b1.fall*meanVal.fall.(predNames{1}) + mdl.fall.(respNames{i}).Coefficients.Estimate(j+1)*x...
                    + b3.fall*meanVal.fall.(predNames{3}) + b4.fall*meanVal.fall.(predNames{4})...
                    + b5.fall*meanVal.fall.(predNames{5}),'k')
            case 'chla'
                plot(x, b0.fall + b1.fall*meanVal.fall.(predNames{1}) + b2.fall*meanVal.fall.(predNames{2})...
                    + mdl.fall.(respNames{i}).Coefficients.Estimate(j+1)*x + b4.fall*meanVal.fall.(predNames{4})...
                    + b5.fall*meanVal.fall.(predNames{5}) + b6.fall*meanVal.fall.(predNames{6}),'k')
            case 'turbidity'
                plot(x, b0.fall + b1.fall*meanVal.fall.(predNames{1}) + b2.fall*meanVal.fall.(predNames{2})...
                    + b3.fall*meanVal.fall.(predNames{3}) + mdl.fall.(respNames{i}).Coefficients.Estimate(j+1)*x...
                    + b5.fall*meanVal.fall.(predNames{5}) + b6.fall*meanVal.fall.(predNames{6}),'k')
            case 'tidal'
                plot(x, b0.fall + b1.fall*meanVal.fall.(predNames{1}) + b2.fall*meanVal.fall.(predNames{2})...
                    + b3.fall*meanVal.fall.(predNames{3}) + b4.fall*meanVal.fall.(predNames{4})...
                    + mdl.fall.(respNames{i}).Coefficients.Estimate(j+1)*x + b6.fall*meanVal.fall.(predNames{6}),'k')
            case 'wspd'
                plot(x, b0.fall + b1.fall*meanVal.fall.(predNames{1}) + b2.fall*meanVal.fall.(predNames{2})...
                    + b3.fall*meanVal.fall.(predNames{3}) + b4.fall*meanVal.fall.(predNames{4})...
                    + b5.fall*meanVal.fall.(predNames{5}) + mdl.fall.(respNames{i}).Coefficients.Estimate(j+1)*x,'k')
        end
        pbaspect(ax4,[1 1 1]); % Make relative lengths of axes equal
        % make a text object for the title
        xl = get(gca(),'Xlim');
        yl = get(gca(),'Ylim');
        text(xl(1),yl(2),'Fall','VerticalAlignment','top','HorizontalAlignment','left','FontSize',12)
        % text(xl(1)+.5,yl(2)-5,'Fall','VerticalAlignment','top','HorizontalAlignment','left','FontSize',12)

        % Control x- and y-tick labels
        ax = t1.Children;
        xlbl = get(ax,'XLabel');
        ylbl = get(ax,'YLabel');
        set([xlbl{:}],'String','');
        set([ylbl{:}],'String','');
        set(ax,'XTickLabel',{})
        set(ax,'YTickLabel',{})
        set([ax(1),ax(2)],'XTickLabelMode','auto')
        set([ax(2),ax(4)],'YTickLabelMode','auto')

        xlabel(t1,xlblNames(j),'FontSize',14)
    end
    ylabel(t,strcat(respNames(i),' (mmol O_2 m^{-2} d^{-1})'),'FontSize',14)
end

%====Save the plots========================================================
option = questdlg('Save plots as .png and .fig?','Save plots','Yes','No','Yes');

switch option
    case 'Yes'
        cd([rootpath,'figures\stats-analyses\metabolism\3x2 Seasonal Plots\All Sites'])
        saveas(fig(1),'GPP.png')
        saveas(fig(1),'GPP.fig')
        saveas(fig(2),'ER.png')
        saveas(fig(2),'ER.fig')
        saveas(fig(3),'NEM.png')
        saveas(fig(3),'NEM.fig')
        disp('Plots saved!')
    case 'No'
        disp('Plots not saved.')
end

%%
% Create a table for ease of copying results to manuscript

mdl = allSites_mdl_norm.summer; % Change the seasons

for i = 1:length(respNames)
    coeffs.(respNames{i}) = table('RowNames',{'Estimate','p-value'});
    coeffs.(respNames{i}).S = [mdl.(respNames{i}).Coefficients.Estimate(2); mdl.(respNames{i}).Coefficients.pValue(2)];
    coeffs.(respNames{i}).T = [mdl.(respNames{i}).Coefficients.Estimate(3); mdl.(respNames{i}).Coefficients.pValue(3)];
    coeffs.(respNames{i}).chla = [mdl.(respNames{i}).Coefficients.Estimate(4); mdl.(respNames{i}).Coefficients.pValue(4)];
    coeffs.(respNames{i}).turbidity = [mdl.(respNames{i}).Coefficients.Estimate(5); mdl.(respNames{i}).Coefficients.pValue(5)];
    coeffs.(respNames{i}).tidal = [mdl.(respNames{i}).Coefficients.Estimate(6); mdl.(respNames{i}).Coefficients.pValue(6)];
    coeffs.(respNames{i}).wspd = [mdl.(respNames{i}).Coefficients.Estimate(7); mdl.(respNames{i}).Coefficients.pValue(7)];
end

    %%
    %====MLRs of site average: (1) for each season and (2) across all seasons==
    dat_syn = synchronize(north_TT,south_TT,gull_TT,'daily');
    mean_depth = mean([dat_syn.depth_north_TT,dat_syn.depth_gull_TT,dat_syn.depth_south_TT],2,'omitmissing');

    % meanSites_dat =

    %%
    %===MLR for Gull by season=================================================

    %----Original (not scaled) data--------------------------------------------
    % Break data into seasons
    gull_dat.winter = gull_dat.all(indWinter,:);
    gull_dat.spring = gull_dat.all(indSpring,:);
    gull_dat.summer = gull_dat.all(indSummer,:);
    gull_dat.fall = gull_dat.all(indFall,:);

    % Compute models with original data
    gull_mdl.winter.GPP = fitlm(gull_dat.winter,'GPP ~ depth + wspd + temperature + salinity + chla + turbidity');
    gull_mdl.spring.GPP = fitlm(gull_dat.spring,'GPP ~ depth + wspd + temperature + salinity + chla + turbidity');
    gull_mdl.summer.GPP = fitlm(gull_dat.summer,'GPP ~ depth + wspd + temperature + salinity + chla + turbidity');
    gull_mdl.fall.GPP = fitlm(gull_dat.fall,'GPP ~ depth + wspd + temperature + salinity + chla + turbidity');

    gull_mdl.winter.ER = fitlm(gull_dat.winter,'ER ~ depth + wspd + temperature + salinity + chla + turbidity');
    gull_mdl.spring.ER = fitlm(gull_dat.spring,'ER ~ depth + wspd + temperature + salinity + chla + turbidity');
    gull_mdl.summer.ER = fitlm(gull_dat.summer,'ER ~ depth + wspd + temperature + salinity + chla + turbidity');
    gull_mdl.fall.ER = fitlm(gull_dat.fall,'ER ~ depth + wspd + temperature + salinity + chla + turbidity');

    gull_mdl.winter.NEM = fitlm(gull_dat.winter,'NEM ~ depth + wspd + temperature + salinity + chla + turbidity');
    gull_mdl.spring.NEM = fitlm(gull_dat.spring,'NEM ~ depth + wspd + temperature + salinity + chla + turbidity');
    gull_mdl.summer.NEM = fitlm(gull_dat.summer,'NEM ~ depth + wspd + temperature + salinity + chla + turbidity');
    gull_mdl.fall.NEM = fitlm(gull_dat.fall,'NEM ~ depth + wspd + temperature + salinity + chla + turbidity');

    %----Scaled data-----------------------------------------------------------
    % Scale response and explanatory variables to center 0 and std 1 (Lowe et al., 2019)
    gull_dat_norm.all = normalize(gull_dat.all(:,1:13)); % Normalize all parameters except "season" and "site"
    gull_dat_norm.all.("season") = gull_dat.all.season; % Add "season" column back in

    % Break normalized data into seasons
    gull_dat_norm.winter = gull_dat_norm.all(indWinter,:);
    gull_dat_norm.spring = gull_dat_norm.all(indSpring,:);
    gull_dat_norm.summer = gull_dat_norm.all(indSummer,:);
    gull_dat_norm.fall = gull_dat_norm.all(indFall,:);

    % Compute models with scaled data
    gull_mdl_norm.winter.GPP = fitlm(gull_dat_norm.winter,'GPP ~ depth + wspd + temperature + salinity + chla + turbidity');
    gull_mdl_norm.spring.GPP = fitlm(gull_dat_norm.spring,'GPP ~ depth + wspd + temperature + salinity + chla + turbidity');
    gull_mdl_norm.summer.GPP = fitlm(gull_dat_norm.summer,'GPP ~ depth + wspd + temperature + salinity + chla + turbidity');
    gull_mdl_norm.fall.GPP = fitlm(gull_dat_norm.fall,'GPP ~ depth + wspd + temperature + salinity + chla + turbidity');

    gull_mdl_norm.winter.ER = fitlm(gull_dat_norm.winter,'ER ~ depth + wspd + temperature + salinity + chla + turbidity');
    gull_mdl_norm.spring.ER = fitlm(gull_dat_norm.spring,'ER ~ depth + wspd + temperature + salinity + chla + turbidity');
    gull_mdl_norm.summer.ER = fitlm(gull_dat_norm.summer,'ER ~ depth + wspd + temperature + salinity + chla + turbidity');
    gull_mdl_norm.fall.ER = fitlm(gull_dat_norm.fall,'ER ~ depth + wspd + temperature + salinity + chla + turbidity');

    gull_mdl_norm.winter.NEM = fitlm(gull_dat_norm.winter,'NEM ~ depth + wspd + temperature + salinity + chla + turbidity');
    gull_mdl_norm.spring.NEM = fitlm(gull_dat_norm.spring,'NEM ~ depth + wspd + temperature + salinity + chla + turbidity');
    gull_mdl_norm.summer.NEM = fitlm(gull_dat_norm.summer,'NEM ~ depth + wspd + temperature + salinity + chla + turbidity');
    gull_mdl_norm.fall.NEM = fitlm(gull_dat_norm.fall,'NEM ~ depth + wspd + temperature + salinity + chla + turbidity');
