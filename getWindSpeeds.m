%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% getWindSpeeds.m
% This script creates output tables (.mat) of wind speed data (1) measured at
% the Gull Island met station and (2) downloaded from ERA5.
%
% AUTHOR:
% Emily Chua 
% 
% DATE:
% First created: 4/2/2024
% Last Updated: 4/3/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;close all;clc

rootpath = 'G:\My Drive\Postdoc\Work\SMIIL\';

cd([rootpath,'physical-data\wind-speed\gull-met'])

%====Import the Gull Met wind speed data===================================
dat1 = readtable('Gull_Met_210914_220625.csv');
varNames = ["datetime_utc","Tair","wspd","wdir","patm","rhumid","source"];
varUnits = ["","degC","m/s","","hPa","",""];
dat1.Properties.VariableNames = varNames;
dat1.Properties.VariableUnits = varUnits;

dat2 = readtable('Gull_Met_WL_220625_220824.xlsx');
varNames = ["datetime_utc","wspd","wdir","Tair","rhumid","patm"];
varUnits = ["","m/s","","degC","","hPa"];
dat2.Properties.VariableNames = varNames;
dat2.Properties.VariableUnits = varUnits;
dat2.datetime_utc = datetime(dat2.datetime_utc);

dat3 = readtable('Gull_Met_WL_220824_221008.xlsx');
dat3.Properties.VariableNames = varNames;
dat3.Properties.VariableUnits = varUnits;
dat3.datetime_utc = datetime(dat3.datetime_utc);

dat4 = readtable('Gull_Met_WL_221008_to_230207.xlsx');
dat4.Properties.VariableNames = varNames;
dat4.Properties.VariableUnits = varUnits;
dat4.datetime_utc = datetime(dat4.datetime_utc);

dat5 = readtable('SMIIL_Gull_Met_WL_230207_to_230524.xlsx');
dat5.Properties.VariableNames = varNames;
dat5.Properties.VariableUnits = varUnits;
dat5.datetime_utc = datetime(dat5.datetime_utc);

% Concatenate data tables
dat1 = dat1(:,{'datetime_utc','wspd','wdir','Tair','rhumid','patm'});   % Restructure table 1
metDat_orig = [dat1;dat2;dat3;dat4;dat5];
metDat_orig.datetime_utc.TimeZone = 'UTC';

metDat_orig = table2timetable(metDat_orig);

clear dat1 dat2 dat3 dat4 dat5 varNames varUnits

%====QC the Gull Met data=======================================================
% Visually assess the data
fig1 = figure(1);clf
fig1.WindowState = 'maximized';
t = tiledlayout(3,1);
ax1 = nexttile;
plot(metDat_orig.datetime_utc,metDat_orig.patm,'.k');ylabel('p_{atm}')
ax2 = nexttile;
plot(metDat_orig.datetime_utc,metDat_orig.Tair,'.k');ylabel('T_{air} (^oC)')
ax3 = nexttile;
plot(metDat_orig.datetime_utc,metDat_orig.wspd,'.k');ylabel('Wind speed (m/s)')
title(t,'Gull Island Met Station - Original Data','fontsize',16)
t.TileSpacing = 'compact';

% Link the axes
linkaxes([ax1,ax2,ax3],'x')     

% Create a table of flags for each sonde
metDat_cleaned = metDat_orig;

% Start with everything as passing (flag = 1)
flags = ones(height(metDat_orig),width(metDat_orig));
flags = array2table(flags);
varNames = ["Tair","wspd","wdir","patm","rhumid"];
flags.Properties.VariableNames = varNames;

% Gross Range Test - Atmospheric Pressure 
% INPUTS
p_low = 900;      % Lower limit (hPa)
p_high = 1100;    % Upper limit (hPa)

ind_low = find(metDat_orig.patm < p_low);
ind_high = find(metDat_orig.patm > p_high);
ind_grossRange = [ind_low;ind_high];

% Gross Range Test - Wind Speed
% INPUTS
u_low = 0;      % Lower limit (m/s)
u_high = 25;    % Upper limit (m/s) -- highest Tropical Storm Ian wind speed was 22.7 m/s

ind_low = find(metDat_orig.wspd < u_low);
ind_high = find(metDat_orig.wspd > u_high);
ind_grossRange = [ind_grossRange;ind_low;ind_high];

% Clean data
% Discard points that failed Gross Range Test
metDat_cleaned(ind_grossRange,:) = {NaN};

% Flag all failed points
flags(ind_grossRange,:) = {4};

cd([rootpath,'figures\diel-analysis-figures\physical-data'])

% Highlight the flagged points
fig2 = figure(2);clf
fig2.WindowState = 'maximized';
t = tiledlayout(3,1);
ax1 = nexttile;
plot(metDat_orig.datetime_utc,metDat_orig.patm,'.k','MarkerSize',6,'DisplayName','Original Data')
ylabel('p_{atm}')
hold on
plot(metDat_orig.datetime_utc(ind_grossRange),metDat_orig.patm(ind_grossRange),'or','MarkerSize',8,'DisplayName','Flagged Points')
legend('show')
ax2 = nexttile;
plot(metDat_orig.datetime_utc,metDat_orig.Tair,'.k','MarkerSize',6);ylabel('T_{air} (^oC)')
hold on
plot(metDat_orig.datetime_utc(ind_grossRange),metDat_orig.Tair(ind_grossRange),'or','MarkerSize',8)
ax3 = nexttile;
plot(metDat_orig.datetime_utc,metDat_orig.wspd,'.k','MarkerSize',6);ylabel('Wind speed (m/s)')
hold on
plot(metDat_orig.datetime_utc(ind_grossRange),metDat_orig.wspd(ind_grossRange),'or','MarkerSize',8)
title(t,'Gull Island Met Station - Flagged Points','fontsize',16)
t.TileSpacing = 'compact';
linkaxes([ax1,ax2,ax3],'x')     

% Plot the cleaned data
fig3 = figure(3);clf
fig3.WindowState = 'maximized';
t = tiledlayout(3,1);
ax1 = nexttile;
plot(metDat_cleaned.datetime_utc,metDat_cleaned.patm,'.k','MarkerSize',6);ylabel('p_{atm}')
hold on
plot(metDat_cleaned.datetime_utc(ind_grossRange),metDat_cleaned.patm(ind_grossRange),'or','MarkerSize',8)
ax2 = nexttile;
plot(metDat_cleaned.datetime_utc,metDat_cleaned.Tair,'.k','MarkerSize',6);ylabel('T_{air} (^oC)')
hold on
plot(metDat_cleaned.datetime_utc(ind_grossRange),metDat_cleaned.Tair(ind_grossRange),'or','MarkerSize',8)
ax3 = nexttile;
plot(metDat_cleaned.datetime_utc,metDat_cleaned.wspd,'.k','MarkerSize',6);ylabel('Wind speed (m/s)')
hold on
plot(metDat_cleaned.datetime_utc(ind_grossRange),metDat_cleaned.wspd(ind_grossRange),'or','MarkerSize',8)
title(t,'Gull Island Met Station - Cleaned Data','fontsize',16)
t.TileSpacing = 'compact';
linkaxes([ax1,ax2,ax3],'x')     

%====Import wind data downloaded from ERA5=================================
% https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels?tab=form
cd([rootpath,'\physical-data\wind-speed\ERA5'])

fileName1 = 'wind_2021.nc';
timeData = ncread(fileName1,'time');
timeData1 = double(timeData);
u10 = ncread(fileName1,'u10');
v10 = ncread(fileName1,'v10');
for i = 1:length(timeData1)
    u10_1(i,1) = u10(1,1,i);
    v10_1(i,1) = v10(1,1,i);
end

fileName2 = 'wind_2022.nc';
timeData = ncread(fileName2,'time');
timeData2 = double(timeData);
u10 = ncread(fileName2,'u10');
v10 = ncread(fileName2,'v10');
for i = 1:length(timeData2)
    u10_2(i,1) = u10(1,1,i);
    v10_2(i,1) = v10(1,1,i);
end

fileName3 = 'wind_2023.nc';
timeData = ncread(fileName3,'time');
timeData3 = double(timeData);
u10 = ncread(fileName3,'u10');
v10 = ncread(fileName3,'v10');
for i = 1:length(timeData3)
    u10_3(i,1) = u10(1,1,i);
    v10_3(i,1) = v10(1,1,i);
end

fileName4 = 'wind_2024.nc';
timeData = ncread(fileName4,'time');
timeData4 = double(timeData);
u10 = ncread(fileName4,'u10');
v10 = ncread(fileName4,'v10');
for i = 1:length(timeData4)
    u10_4(i,1) = u10(1,1,1,i);
    v10_4(i,1) = v10(1,1,1,i);
end
lat = ncread(fileName4,'latitude');
lon = ncread(fileName4,'longitude');

clearvars timeData u10 v10

u10 = [u10_1;u10_2;u10_3;u10_4];
v10 = [v10_1;v10_2;v10_3;v10_4];
timeData = [timeData1;timeData2;timeData3;timeData4];

%====Convert time unit=====================================================
% (https://memg.ocean.dal.ca/grosse/Teaching/MM2019/Materials/Labs/Lab10/Lab_10.pdf)
timeUnit = ncreadatt(fileName1,'time','units');

switch timeUnit(1:3)
    case 'sec'
        timeFac = 1/86400;
    case 'min'
        timeFac = 1/1440;
    case {'hou', 'hrs'}
        timeFac = 1/24;
    case 'day'
        timeFac = 1;
    otherwise
        error('Invalid time unit.');
end

% get the reference time used in the NetCDF file
i = strfind(timeUnit, 'since') + length('since');
refTimeStr = strtrim(timeUnit(i:end));
refTimeNum = datenum(refTimeStr);

% do time conversion
timeData = refTimeNum + timeData*timeFac;

% convert the day counter into a date vector
timeVec = datevec(timeData);
datetime = datetime(timeVec);
datetime.TimeZone = 'UTC';

%====Calculate wind speed from u- and v-components=========================
% (https://confluence.ecmwf.int/pages/viewpage.action?pageId=133262398)
wspd = sqrt(u10.^2 + v10.^2);

%====Create time table of ERA5 data========================================
era5Dat = timetable(datetime,wspd);
era5Dat.Properties.VariableUnits = {'m/s'};

%====Compare cleaned Gull met data and ERA5 data===========================
cd([rootpath,'figures\diel-analysis-figures\physical-data'])

fig4 = figure(4);clf
fig4.WindowState = 'maximized';
plot(metDat_cleaned.datetime_utc,metDat_cleaned.wspd,'.k','DisplayName','Gull Met Station')
hold on
plot(era5Dat.datetime,era5Dat.wspd,'.','DisplayName','ERA5')
hold off
ylabel('Wind speed (m/s)')
legend('show')

% Compare Gull and ERA5 hourly mean wind speeds
metDat_hourlyMean = retime(metDat_cleaned,'hourly','mean');

% Find data range in ERA5 data that spans Gull Met data
[~,start] = ismember(metDat_hourlyMean.datetime_utc(1),era5Dat.datetime);
[~,stop] = ismember(metDat_hourlyMean.datetime_utc(end),era5Dat.datetime);

fig5 = figure(5);clf
fig5.WindowState = 'maximized';
plot(metDat_hourlyMean.datetime_utc,metDat_hourlyMean.wspd,'.k','DisplayName','Gull Met Station - Hourly Means')
hold on
plot(era5Dat.datetime,era5Dat.wspd,'.','DisplayName','ERA5')
hold off
ylabel('Wind speed (m/s)')
legend('show')

% Plot linear regression between cleaned Gull Met and ERA5 data
tbl = table(metDat_hourlyMean.wspd,era5Dat.wspd(start:stop));
tbl.Properties.VariableNames = ["metDat","era5"];
mdl = fitlm(tbl.metDat,tbl.era5,'y~x1-1');  % Force intercept through zero; see Wilkinson notation

% Create equation string for plot
eqn = ['y = ',num2str(mdl.Coefficients.Estimate,3),'x'];
R2 = num2str(mdl.Rsquared.Ordinary,2);

fig6 = figure(6);clf
% fig6.WindowState = 'maximized';
h = plot(mdl,'marker','.');
delete(h([3 4]))    % Delete confidence bounds on plot
hold on
plot([0 ylim], [0 ylim],'--k')
legend('Data',[eqn,newline,'R^2 = ',R2],'1:1 line','Location','southeast')
xlabel('Gull Met Station (Cleaned)')
ylabel('ERA5')
title('Hourly Wind Speed (m/s)')
daspect([1 1 1])

%%
%====Save the cleaned Gull Met data and ERA5 data==========================
cd('G:\My Drive\Postdoc\Work\SMIIL\physical-data\wind-speed')
save('windSpeed.mat','metDat_cleaned','era5Dat')
