%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotTimeseries.m
% This script...
%
% AUTHOR:
% Emily Chua
%
% DATE:
% First created: 8/6/2024
% Last updated: 8/11/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;close all;clc

rootpath = 'G:\My Drive\Postdoc\Work\SMIIL\';

%==========================================================================
% Import data
%==========================================================================
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

%====Create tables containing all data for each site=======================
% Gull
dt2 = dateshift(gull_metab.daystart_dt,'start','day');
gull_metab.daystart_dt = dt2;
daily_gull = synchronize(gull_dailyAvg,wspd_dailyAvg,par_dailyAvg,gull_metab);
daily_gull = removevars(daily_gull,{'deployment','dayend_dt','daylength','R_hourly','P_hourly','R_daily','P_daily','datetime_local','Tair','light_lux'});
daily_gull(daily_gull.datetime_utc < datetime(2021,6,29,"TimeZone","UTC"),:) = [];
daily_gull(daily_gull.datetime_utc > datetime(2024,6,4,"TimeZone","UTC"),:) = [];

% North
dt2 = dateshift(north_metab.daystart_dt,'start','day');
north_metab.daystart_dt = dt2;
daily_north = synchronize(north_dailyAvg,wspd_dailyAvg,par_dailyAvg,north_metab);
daily_north = removevars(daily_north,{'deployment','dayend_dt','daylength','R_hourly','P_hourly','R_daily','P_daily','datetime_local','Tair','light_lux'});
daily_north(daily_north.datetime_utc < datetime(2021,6,29,"TimeZone","UTC"),:) = [];
daily_north(daily_north.datetime_utc > datetime(2024,6,4,"TimeZone","UTC"),:) = [];

% South
dt2 = dateshift(south_metab.daystart_dt,'start','day');
south_metab.daystart_dt = dt2;
daily_south = synchronize(south_dailyAvg,wspd_dailyAvg,par_dailyAvg,south_metab);
daily_south = removevars(daily_south,{'deployment','dayend_dt','daylength','R_hourly','P_hourly','R_daily','P_daily','datetime_local','Tair','light_lux'});
daily_south(daily_south.datetime_utc < datetime(2021,6,29,"TimeZone","UTC"),:) = [];
daily_south(daily_south.datetime_utc > datetime(2024,6,4,"TimeZone","UTC"),:) = [];

%%
%==========================================================================
% Perform ANOVA
%==========================================================================

% Create a grouping variable
group = {'North','Gull','South'};

%====Perform one-way ANOVA for GPP=========================================
data = [daily_north.GPP, daily_gull.GPP, daily_south.GPP];
[p, tbl, stats] = anova1(data, group);

% Perform multiple comparisons if the ANOVA is significant
if p < 0.05
    [c,m,h,gnames] = multcompare(stats);
    tbl1 = array2table(c,"VariableNames",["Group A","Group B","Lower Limit","A-B","Upper Limit","P-values"]);
    tbl1.("Group A") = gnames(tbl1.("Group A"));
    tbl1.("Group B") = gnames(tbl1.("Group B"));
end

%====Perform one-way ANOVA for ER=========================================
data = [daily_north.ER, daily_gull.ER, daily_south.ER];
[p, tbl, stats] = anova1(data, group);

% Perform multiple comparisons if the ANOVA is significant
if p < 0.05
    [c,m,h,gnames] = multcompare(stats);
    tbl2 = array2table(c,"VariableNames",["Group A","Group B","Lower Limit","A-B","Upper Limit","P-values"]);
    tbl2.("Group A") = gnames(tbl2.("Group A"));
    tbl2.("Group B") = gnames(tbl2.("Group B"));
end

%====Perform one-way ANOVA for NEM=========================================
data = [daily_north.NEM, daily_gull.NEM, daily_south.NEM];
[p, tbl, stats] = anova1(data, group);

% Perform multiple comparisons if the ANOVA is significant
if p < 0.05
    [c,m,h,gnames] = multcompare(stats);
    tbl3 = array2table(c,"VariableNames",["Group A","Group B","Lower Limit","A-B","Upper Limit","P-values"]);
    tbl3.("Group A") = gnames(tbl3.("Group A"));
    tbl3.("Group B") = gnames(tbl3.("Group B"));
end

%% Old
% Environmental drivers
% data = [daily_north.temperature, daily_gull.temperature, daily_south.temperature];
% data = [daily_north.salinity, daily_gull.salinity, daily_south.salinity];
% data = [daily_north.pH, daily_gull.pH, daily_south.pH];
% data = [daily_north.DOsat, daily_gull.DOsat, daily_south.DOsat];
% data = [daily_north.chla, daily_gull.chla, daily_south.chla];
% data = [daily_north.turbidity, daily_gull.turbidity, daily_south.turbidity];