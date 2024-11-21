%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% conductAnova.m
% This script performs (1) One-way ANOVAs to see if mean GPP, ER, and NEM are
% different depending on site, and (2) Pairwise comparisons using a
% multiple comparison test to identify the sites that have significantly
% different means.
%
% AUTHOR:
% Emily Chua
%
% DATE:
% First created: 8/6/2024
% Last updated: 8/19/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;close all;clc

rootpath = 'G:\My Drive\Postdoc\Work\SMIIL\';

%==========================================================================
% Import data
%==========================================================================
site = 'gull';
cd([rootpath,'open-water-platform-data\',site,'\cleaned\final-qc'])
load([site,'-cleaned.mat']);
params.gull = finalQC;
% cd([rootpath,'diel-method\matlab-results\final-qc\',site])
% load('diel_res.mat')
% metab.gull = diel_dtd;
cd([rootpath,'diel-method\uncertainty-analysis\',site])
load('MonteCarloResults.mat')
metab.gull = diel_dtd_MC;
% Remove times with anomalous GPP and ER values
anomER = find(metab.gull.ER_avg > 0);
anomGPP = find(metab.gull.GPP_avg < 0);
metab.gull([anomER;anomGPP],:) = [];

site = 'north';
cd([rootpath,'open-water-platform-data\',site,'\cleaned\final-qc'])
load([site,'-cleaned.mat']);
params.north = finalQC;
% cd([rootpath,'diel-method\matlab-results\final-qc\',site])
% load('diel_res.mat')
% metab.north = diel_dtd;
cd([rootpath,'diel-method\uncertainty-analysis\',site])
load('MonteCarloResults.mat')
metab.north = diel_dtd_MC;
% Remove times with anomalous GPP and ER values
anomER = find(metab.north.ER_avg > 0);
anomGPP = find(metab.north.GPP_avg < 0);
metab.north([anomER;anomGPP],:) = [];

site = 'south';
cd([rootpath,'open-water-platform-data\',site,'\cleaned\final-qc'])
load([site,'-cleaned.mat']);
params.south = finalQC;
% cd([rootpath,'diel-method\matlab-results\final-qc\',site])
% load('diel_res.mat')
% metab.south = diel_dtd;
cd([rootpath,'diel-method\uncertainty-analysis\',site])
load('MonteCarloResults.mat')
metab.south= diel_dtd_MC;
% Remove times with anomalous GPP and ER values
anomER = find(metab.south.ER_avg > 0);
anomGPP = find(metab.south.GPP_avg < 0);
metab.south([anomER;anomGPP],:) = [];

clearvars finalQC diel_dtd diel_obs

%====Import windspeed and PAR data=========================================
cd([rootpath,'physical-data\final-dataset'])
load('windspeed.mat')

cd([rootpath,'physical-data\final-dataset'])
load('par.mat')

%====Calculate daily means for each site===================================
gull_dailyAvg = retime(params.gull,'daily','mean');
north_dailyAvg = retime(params.north,'daily','mean');
south_dailyAvg = retime(params.south,'daily','mean');

wspd_dailyAvg = retime(era5Dat,'daily','mean');
par_dailyAvg= retime(parDat,'daily','mean');

cd([rootpath,'figures\stats-analyses'])

%====Create tables containing all data for each site=======================
% Gull
dt2 = dateshift(metab.gull.daystart_dt,'start','day');
metab.gull.daystart_dt = dt2;
daily_gull = synchronize(gull_dailyAvg,wspd_dailyAvg,par_dailyAvg,metab.gull);
daily_gull = removevars(daily_gull,{'deployment','dayend_dt','daylength','R_hourly_avg','P_hourly_avg','R_daily_avg','P_daily_avg','datetime_local','Tair','light_lux'});
daily_gull(daily_gull.datetime_utc < datetime(2021,6,29,"TimeZone","UTC"),:) = [];
daily_gull(daily_gull.datetime_utc > datetime(2024,6,4,"TimeZone","UTC"),:) = [];

% North
dt2 = dateshift(metab.north.daystart_dt,'start','day');
metab.north.daystart_dt = dt2;
daily_north = synchronize(north_dailyAvg,wspd_dailyAvg,par_dailyAvg,metab.north);
daily_north = removevars(daily_north,{'deployment','dayend_dt','daylength','R_hourly_avg','P_hourly_avg','R_daily_avg','P_daily_avg','datetime_local','Tair','light_lux'});
daily_north(daily_north.datetime_utc < datetime(2021,6,29,"TimeZone","UTC"),:) = [];
daily_north(daily_north.datetime_utc > datetime(2024,6,4,"TimeZone","UTC"),:) = [];

% South
dt2 = dateshift(metab.south.daystart_dt,'start','day');
metab.south.daystart_dt = dt2;
daily_south = synchronize(south_dailyAvg,wspd_dailyAvg,par_dailyAvg,metab.south);
daily_south = removevars(daily_south,{'deployment','dayend_dt','daylength','R_hourly_avg','P_hourly_avg','R_daily_avg','P_daily_avg','datetime_local','Tair','light_lux'});
daily_south(daily_south.datetime_utc < datetime(2021,6,29,"TimeZone","UTC"),:) = [];
daily_south(daily_south.datetime_utc > datetime(2024,6,4,"TimeZone","UTC"),:) = [];

%==========================================================================
% Perform ANOVA
%==========================================================================

% Create a grouping variable
group = {'North','Gull','South'};

%====Perform one-way ANOVA for GPP=========================================
data = [daily_north.GPP_avg, daily_gull.GPP_avg, daily_south.GPP_avg];
[p, tbl, stats] = anova1(data, group);

% Perform multiple comparisons if the ANOVA is significant
if p < 0.05
    [c,m,h,gnames] = multcompare(stats);
    tbl1 = array2table(c,"VariableNames",["Group A","Group B","Lower Limit","A-B","Upper Limit","P-values"]);
    tbl1.("Group A") = gnames(tbl1.("Group A"));
    tbl1.("Group B") = gnames(tbl1.("Group B"));
end

%====Perform one-way ANOVA for ER=========================================
data = [daily_north.ER_avg, daily_gull.ER_avg, daily_south.ER_avg];
[p, tbl, stats] = anova1(data, group);

% Perform multiple comparisons if the ANOVA is significant
if p < 0.05
    [c,m,h,gnames] = multcompare(stats);
    tbl2 = array2table(c,"VariableNames",["Group A","Group B","Lower Limit","A-B","Upper Limit","P-values"]);
    tbl2.("Group A") = gnames(tbl2.("Group A"));
    tbl2.("Group B") = gnames(tbl2.("Group B"));
end

%====Perform one-way ANOVA for NEM=========================================
data = [daily_north.NEM_avg, daily_gull.NEM_avg, daily_south.NEM_avg];
[p, tbl, stats] = anova1(data, group);

% Perform multiple comparisons if the ANOVA is significant
if p < 0.05
    [c,m,h,gnames] = multcompare(stats);
    tbl3 = array2table(c,"VariableNames",["Group A","Group B","Lower Limit","A-B","Upper Limit","P-values"]);
    tbl3.("Group A") = gnames(tbl3.("Group A"));
    tbl3.("Group B") = gnames(tbl3.("Group B"));
end
