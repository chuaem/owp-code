%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% statsAnalysis_metabolicDrivers.m
% This script...
%
% AUTHOR:
% Emily Chua
%
% DATE:
% First created: 7/22/2024
% Last updated:
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

%====Calculate monthly means for each site=================================
% Gull
dt2 = dateshift(gull_metab.daystart_dt,'start','day');
gull_metab.daystart_dt = dt2;
TT_gull = synchronize(gull_dailyAvg,wspd_dailyAvg,par_dailyAvg,gull_metab);
TT_gull = removevars(TT_gull,{'deployment','dayend_dt','daylength','R_hourly','P_hourly','R_daily','P_daily','datetime_local','Tair','light_lux'});
TT_gull(TT_gull.datetime_utc < datetime(2021,6,29,"TimeZone","UTC"),:) = [];
TT_gull(TT_gull.datetime_utc > datetime(2024,6,4,"TimeZone","UTC"),:) = [];

% Get month numbers
mo = month(TT_gull.datetime_utc);

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
TT_gull.season(indWinter) = 1;
TT_gull.season(indSpring) = 2;
TT_gull.season(indSummer) = 3;
TT_gull.season(indFall) = 4;
TT_gull.season = categorical(TT_gull.season);

% Add a column for site number
TT_gull.site = repelem(2,height(TT_gull))';
TT_gull.site = categorical(TT_gull.site);

% North
dt2 = dateshift(north_metab.daystart_dt,'start','day');
north_metab.daystart_dt = dt2;
TT_north = synchronize(north_dailyAvg,wspd_dailyAvg,par_dailyAvg,north_metab);
TT_north = removevars(TT_north,{'deployment','dayend_dt','daylength','R_hourly','P_hourly','R_daily','P_daily','datetime_local','Tair','light_lux'});
TT_north(TT_north.datetime_utc < datetime(2021,6,29,"TimeZone","UTC"),:) = [];
TT_north(TT_north.datetime_utc > datetime(2024,6,4,"TimeZone","UTC"),:) = [];

% Get month numbers
mo = month(TT_north.datetime_utc);

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
TT_north.season(indWinter) = 1;
TT_north.season(indSpring) = 2;
TT_north.season(indSummer) = 3;
TT_north.season(indFall) = 4;
TT_north.season = categorical(TT_north.season);

% Add a column for site number
TT_north.site = repelem(1,height(TT_north))';
TT_north.site = categorical(TT_north.site);

% South
dt2 = dateshift(south_metab.daystart_dt,'start','day');
south_metab.daystart_dt = dt2;
TT_south = synchronize(south_dailyAvg,wspd_dailyAvg,par_dailyAvg,south_metab);
TT_south = removevars(TT_south,{'deployment','dayend_dt','daylength','R_hourly','P_hourly','R_daily','P_daily','datetime_local','Tair','light_lux'});
TT_south(TT_south.datetime_utc < datetime(2021,6,29,"TimeZone","UTC"),:) = [];
TT_south(TT_south.datetime_utc > datetime(2024,6,4,"TimeZone","UTC"),:) = [];

% Get month numbers
mo = month(TT_south.datetime_utc);

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
TT_south.season(indWinter) = 1;
TT_south.season(indSpring) = 2;
TT_south.season(indSummer) = 3;
TT_south.season(indFall) = 4;
TT_south.season = categorical(TT_south.season);

% Add a column for site number
TT_south.site = repelem(3,height(TT_south))';
TT_south.site = categorical(TT_south.site);

%====Create a table to store predictor variables for all sites=============
gull_dat.all = timetable2table(TT_gull);
gull_dat.all.datetime_utc = [];

north_dat.all = timetable2table(TT_north);
north_dat.all.datetime_utc = [];

south_dat.all = timetable2table(TT_south);
south_dat.all.datetime_utc = [];

allSites_tbl = [gull_dat.all;north_dat.all;south_dat.all];

%==========================================================================
% Construct multiple linear regression models for ER, GPP, and NEM
%==========================================================================

%====MLR for Gull by season================================================
% % Original data
% % Break data into seasons
% gull_dat.winter = gull_dat.all(gull_dat.all.season == '1',:);
% gull_dat.spring = gull_dat.all(gull_dat.all.season == '2',:);
% gull_dat.summer = gull_dat.all(gull_dat.all.season == '3',:);
% gull_dat.fall = gull_dat.all(gull_dat.all.season == '4',:);
% 
% % Compute models with original data
% gull_mdl.winter.GPP = fitlm(gull_dat.winter,'GPP ~ depth + temperature + salinity + chla + turbidity + wspd');
% gull_mdl.spring.GPP = fitlm(gull_dat.spring,'GPP ~ depth + temperature + salinity + chla + turbidity + wspd');
% gull_mdl.summer.GPP = fitlm(gull_dat.summer,'GPP ~ depth + temperature + salinity + chla + turbidity + wspd');
% gull_mdl.fall.GPP = fitlm(gull_dat.fall,'GPP ~ depth + temperature + salinity + chla + turbidity + wspd');
% 
% gull_mdl.winter.ER = fitlm(gull_dat.winter,'ER ~ depth + temperature + salinity + chla + turbidity + wspd');
% gull_mdl.spring.ER = fitlm(gull_dat.spring,'ER ~ depth + temperature + salinity + chla + turbidity + wspd');
% gull_mdl.summer.ER = fitlm(gull_dat.summer,'ER ~ depth + temperature + salinity + chla + turbidity + wspd');
% gull_mdl.fall.ER = fitlm(gull_dat.fall,'ER ~ depth + temperature + salinity + chla + turbidity + wspd');
% 
% gull_mdl.winter.NEM = fitlm(gull_dat.winter,'NEM ~ depth + temperature + salinity + chla + turbidity + wspd');
% gull_mdl.spring.NEM = fitlm(gull_dat.spring,'NEM ~ depth + temperature + salinity + chla + turbidity + wspd');
% gull_mdl.summer.NEM = fitlm(gull_dat.summer,'NEM ~ depth + temperature + salinity + chla + turbidity + wspd');
% gull_mdl.fall.NEM = fitlm(gull_dat.fall,'NEM ~ depth + temperature + salinity + chla + turbidity + wspd');

% Scaled data
% Scale response and explanatory variables to center 0 and std 1 (Lowe et al., 2019)
gull_dat_norm.all = normalize(gull_dat.all(:,1:13));
gull_dat_norm.all.("season") = gull_dat.all.season;
% Break data into seasons
gull_dat_norm.winter = gull_dat_norm.all(gull_dat_norm.all.season == '1',:);
gull_dat_norm.spring = gull_dat_norm.all(gull_dat_norm.all.season == '2',:);
gull_dat_norm.summer = gull_dat_norm.all(gull_dat_norm.all.season == '3',:);
gull_dat_norm.fall = gull_dat_norm.all(gull_dat_norm.all.season == '4',:);

% Compute models with scaled data
gull_mdl_norm.winter.GPP = fitlm(gull_dat_norm.winter,'GPP ~ depth + temperature + salinity + chla + turbidity + wspd');
gull_mdl_norm.spring.GPP = fitlm(gull_dat_norm.spring,'GPP ~ depth + temperature + salinity + chla + turbidity + wspd');
gull_mdl_norm.summer.GPP = fitlm(gull_dat_norm.summer,'GPP ~ depth + temperature + salinity + chla + turbidity + wspd');
gull_mdl_norm.fall.GPP = fitlm(gull_dat_norm.fall,'GPP ~ depth + temperature + salinity + chla + turbidity + wspd');

gull_mdl_norm.winter.ER = fitlm(gull_dat_norm.winter,'ER ~ depth + temperature + salinity + chla + turbidity + wspd');
gull_mdl_norm.spring.ER = fitlm(gull_dat_norm.spring,'ER ~ depth + temperature + salinity + chla + turbidity + wspd');
gull_mdl_norm.summer.ER = fitlm(gull_dat_norm.summer,'ER ~ depth + temperature + salinity + chla + turbidity + wspd');
gull_mdl_norm.fall.ER = fitlm(gull_dat_norm.fall,'ER ~ depth + temperature + salinity + chla + turbidity + wspd');

gull_mdl_norm.winter.NEM = fitlm(gull_dat_norm.winter,'NEM ~ depth + temperature + salinity + chla + turbidity + wspd');
gull_mdl_norm.spring.NEM = fitlm(gull_dat_norm.spring,'NEM ~ depth + temperature + salinity + chla + turbidity + wspd');
gull_mdl_norm.summer.NEM = fitlm(gull_dat_norm.summer,'NEM ~ depth + temperature + salinity + chla + turbidity + wspd');
gull_mdl_norm.fall.NEM = fitlm(gull_dat_norm.fall,'NEM ~ depth + temperature + salinity + chla + turbidity + wspd');



%====MLR for all data across all 4 seasons=================================
% mdl.ER = fitlm(allSites_tbl,'ER ~ temperature + turbidity + chla + wspd + season + site');
%
% mdl.GPP = fitlm(allSites_tbl,'GPP ~ temperature + turbidity + chla + wspd + season + site');
%
% mdl.NEM = fitlm(allSites_tbl,'NEM ~ temperature + turbidity + chla + wspd + season + site');

%%
%==========================================================================
% Plot results
%==========================================================================
data = gull_dat;    % Choose which site to plot

%====ER vs. Explanatory Variables==========================================
fig = figure(1);clf
fig.Position = [363.4 905.0 676.8 957.6];
t = tiledlayout(3,2,'Padding','none','TileSpacing','loose');

%----Tile 1: Depth---------------------------------------------------------
t1 = tiledlayout(t,2,2,'Tilespacing','none');
t1.Layout.Tile = 1;

ax1 = nexttile(t1);
plot(data.all.depth,data.all.ER,'.','color',rgb('lightgrey'))
hold on
plot(data.winter.depth,data.winter.ER,'.k')
ylim([-300 100])
pbaspect(ax1,[1 1 1]); % Make relative lengths of axes equal
% make a text object for the title
xl = get(gca(),'Xlim');
yl = get(gca(),'Ylim');
text(xl(1),yl(2),'Winter','VerticalAlignment','top','HorizontalAlignment','left','FontSize',12)

ax2 = nexttile(t1);
plot(data.all.depth,data.all.ER,'.','color',rgb('lightgrey'))
hold on
plot(data.spring.depth,data.spring.ER,'.k')
ylim([-300 100])
pbaspect(ax2,[1 1 1]); % Make relative lengths of axes equal
% make a text object for the title
xl = get(gca(),'Xlim');
yl = get(gca(),'Ylim');
text(xl(1),yl(2),'Spring','VerticalAlignment','top','HorizontalAlignment','left','FontSize',12)

ax3 = nexttile(t1);
plot(data.all.depth,data.all.ER,'.','color',rgb('lightgrey'))
hold on
plot(data.summer.depth,data.summer.ER,'.k')
ylim([-300 100])
pbaspect(ax3,[1 1 1]); % Make relative lengths of axes equal
% make a text object for the title
xl = get(gca(),'Xlim');
yl = get(gca(),'Ylim');
text(xl(1),yl(2),'Summer','VerticalAlignment','top','HorizontalAlignment','left','FontSize',12)

ax4 = nexttile(t1);
plot(data.all.depth,data.all.ER,'.','color',rgb('lightgrey'))
hold on
plot(data.fall.depth,data.fall.ER,'.k')
ylim([-300 100])
pbaspect(ax4,[1 1 1]); % Make relative lengths of axes equal
% make a text object for the title
xl = get(gca(),'Xlim');
yl = get(gca(),'Ylim');
text(xl(1),yl(2),'Fall','VerticalAlignment','top','HorizontalAlignment','left','FontSize',12)

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

xlabel(t1,'Depth (m)')

%----Tile 2: Wind speed----------------------------------------------------
t2 = tiledlayout(t,2,2,'Tilespacing','none');
t2.Layout.Tile = 2;

ax1 = nexttile(t2);
plot(data.all.wspd,data.all.ER,'.','color',rgb('lightgrey'))
hold on
plot(data.winter.wspd,data.winter.ER,'.k')
ylim([-300 100])
pbaspect(ax1,[1 1 1]); % Make relative lengths of axes equal
% make a text object for the title
xl = get(gca(),'Xlim');
yl = get(gca(),'Ylim');
text(xl(1),yl(2),'Winter','VerticalAlignment','top','HorizontalAlignment','left','FontSize',12)

ax2 = nexttile(t2);
plot(data.all.wspd,data.all.ER,'.','color',rgb('lightgrey'))
hold on
plot(data.spring.wspd,data.spring.ER,'.k')
ylim([-300 100])
pbaspect(ax2,[1 1 1]); % Make relative lengths of axes equal
% make a text object for the title
xl = get(gca(),'Xlim');
yl = get(gca(),'Ylim');
text(xl(1),yl(2),'Spring','VerticalAlignment','top','HorizontalAlignment','left','FontSize',12)

ax3 = nexttile(t2);
plot(data.all.wspd,data.all.ER,'.','color',rgb('lightgrey'))
hold on
plot(data.summer.wspd,data.summer.ER,'.k')
ylim([-300 100])
pbaspect(ax3,[1 1 1]); % Make relative lengths of axes equal
% make a text object for the title
xl = get(gca(),'Xlim');
yl = get(gca(),'Ylim');
text(xl(1),yl(2),'Summer','VerticalAlignment','top','HorizontalAlignment','left','FontSize',12)

ax4 = nexttile(t2);
plot(data.all.wspd,data.all.ER,'.','color',rgb('lightgrey'))
hold on
plot(data.fall.wspd,data.fall.ER,'.k')
ylim([-300 100])
pbaspect(ax4,[1 1 1]); % Make relative lengths of axes equal
% make a text object for the title
xl = get(gca(),'Xlim');
yl = get(gca(),'Ylim');
text(xl(1),yl(2),'Fall','VerticalAlignment','top','HorizontalAlignment','left','FontSize',12)

% Control x- and y-tick labels
ax = t2.Children;
xlbl = get(ax,'XLabel');
ylbl = get(ax,'YLabel');
set([xlbl{:}],'String','');
set([ylbl{:}],'String','');
set(ax,'XTickLabel',{})
set(ax,'YTickLabel',{})
set([ax(1),ax(2)],'XTickLabelMode','auto')
set([ax(2),ax(4)],'YTickLabelMode','auto')

xlabel(t2,'Wind speed (m/s)')

%----Tile 3: Temperature---------------------------------------------------
t3 = tiledlayout(t,2,2,'Tilespacing','none');
t3.Layout.Tile = 3;

ax1 = nexttile(t3);
plot(data.all.temperature,data.all.ER,'.','color',rgb('lightgrey'))
hold on
plot(data.winter.temperature,data.winter.ER,'.k')
ylim([-300 100])
pbaspect(ax1,[1 1 1]); % Make relative lengths of axes equal
% make a text object for the title
xl = get(gca(),'Xlim');
yl = get(gca(),'Ylim');
text(xl(1),yl(2),'Winter','VerticalAlignment','top','HorizontalAlignment','left','FontSize',12)

ax2 = nexttile(t3);
plot(data.all.temperature,data.all.ER,'.','color',rgb('lightgrey'))
hold on
plot(data.spring.temperature,data.spring.ER,'.k')
plot(data.all.temperature,gull_mdl.spring.ER.Coefficients.Estimate(1) + gull_mdl.spring.ER.Coefficients.Estimate(4)*data.all.temperature,'k')
ylim([-300 100])
pbaspect(ax2,[1 1 1]); % Make relative lengths of axes equal
% make a text object for the title
xl = get(gca(),'Xlim');
yl = get(gca(),'Ylim');
text(xl(1),yl(2),'Spring','VerticalAlignment','top','HorizontalAlignment','left','FontSize',12)

ax3 = nexttile(t3);
plot(data.all.temperature,data.all.ER,'.','color',rgb('lightgrey'))
hold on
plot(data.summer.temperature,data.summer.ER,'.k')
ylim([-300 100])
pbaspect(ax3,[1 1 1]); % Make relative lengths of axes equal
% make a text object for the title
xl = get(gca(),'Xlim');
yl = get(gca(),'Ylim');
text(xl(1),yl(2),'Summer','VerticalAlignment','top','HorizontalAlignment','left','FontSize',12)

ax4 = nexttile(t3);
plot(data.all.temperature,data.all.ER,'.','color',rgb('lightgrey'))
hold on
plot(data.fall.temperature,data.fall.ER,'.k')
ylim([-300 100])
pbaspect(ax4,[1 1 1]); % Make relative lengths of axes equal
% make a text object for the title
xl = get(gca(),'Xlim');
yl = get(gca(),'Ylim');
text(xl(1),yl(2),'Fall','VerticalAlignment','top','HorizontalAlignment','left','FontSize',12)

% Control x- and y-tick labels
ax = t3.Children;
xlbl = get(ax,'XLabel');
ylbl = get(ax,'YLabel');
set([xlbl{:}],'String','');
set([ylbl{:}],'String','');
set(ax,'XTickLabel',{})
set(ax,'YTickLabel',{})
set([ax(1),ax(2)],'XTickLabelMode','auto')
set([ax(2),ax(4)],'YTickLabelMode','auto')

xlabel(t3,'Temperature (^oC)')

%----Tile 4: Salinity------------------------------------------------------
t4 = tiledlayout(t,2,2,'Tilespacing','none');
t4.Layout.Tile = 4;

ax1 = nexttile(t4);
plot(data.all.salinity,data.all.ER,'.','color',rgb('lightgrey'))
hold on
plot(data.winter.salinity,data.winter.ER,'.k')
ylim([-300 100])
pbaspect(ax1,[1 1 1]); % Make relative lengths of axes equal
% make a text object for the title
xl = get(gca(),'Xlim');
yl = get(gca(),'Ylim');
text(xl(1),yl(2),'Winter','VerticalAlignment','top','HorizontalAlignment','left','FontSize',12)

ax2 = nexttile(t4);
plot(data.all.salinity,data.all.ER,'.','color',rgb('lightgrey'))
hold on
plot(data.spring.salinity,data.spring.ER,'.k')
ylim([-300 100])
pbaspect(ax2,[1 1 1]); % Make relative lengths of axes equal
% make a text object for the title
xl = get(gca(),'Xlim');
yl = get(gca(),'Ylim');
text(xl(1),yl(2),'Spring','VerticalAlignment','top','HorizontalAlignment','left','FontSize',12)

ax3 = nexttile(t4);
plot(data.all.salinity,data.all.ER,'.','color',rgb('lightgrey'))
hold on
plot(data.summer.salinity,data.summer.ER,'.k')
ylim([-300 100])
pbaspect(ax3,[1 1 1]); % Make relative lengths of axes equal
% make a text object for the title
xl = get(gca(),'Xlim');
yl = get(gca(),'Ylim');
text(xl(1),yl(2),'Summer','VerticalAlignment','top','HorizontalAlignment','left','FontSize',12)

ax4 = nexttile(t4);
plot(data.all.salinity,data.all.ER,'.','color',rgb('lightgrey'))
hold on
plot(data.fall.salinity,data.fall.ER,'.k')
ylim([-300 100])
pbaspect(ax4,[1 1 1]); % Make relative lengths of axes equal
% make a text object for the title
xl = get(gca(),'Xlim');
yl = get(gca(),'Ylim');
text(xl(1),yl(2),'Fall','VerticalAlignment','top','HorizontalAlignment','left','FontSize',12)

% Control x- and y-tick labels
ax = t4.Children;
xlbl = get(ax,'XLabel');
ylbl = get(ax,'YLabel');
set([xlbl{:}],'String','');
set([ylbl{:}],'String','');
set(ax,'XTickLabel',{})
set(ax,'YTickLabel',{})
set([ax(1),ax(2)],'XTickLabelMode','auto')
set([ax(2),ax(4)],'YTickLabelMode','auto')

xlabel(t4,'Salinity')

%----Tile 5: Chl a---------------------------------------------------------
t5 = tiledlayout(t,2,2,'Tilespacing','none');
t5.Layout.Tile = 5;

ax1 = nexttile(t5);
plot(data.all.chla,data.all.ER,'.','color',rgb('lightgrey'))
hold on
plot(data.winter.chla,data.winter.ER,'.k')
ylim([-300 100])
pbaspect(ax1,[1 1 1]); % Make relative lengths of axes equal
% make a text object for the title
xl = get(gca(),'Xlim');
yl = get(gca(),'Ylim');
text(xl(1),yl(2),'Winter','VerticalAlignment','top','HorizontalAlignment','left','FontSize',12)

ax2 = nexttile(t5);
plot(data.all.chla,data.all.ER,'.','color',rgb('lightgrey'))
hold on
plot(data.spring.chla,data.spring.ER,'.k')
ylim([-300 100])
pbaspect(ax2,[1 1 1]); % Make relative lengths of axes equal
% make a text object for the title
xl = get(gca(),'Xlim');
yl = get(gca(),'Ylim');
text(xl(1),yl(2),'Spring','VerticalAlignment','top','HorizontalAlignment','left','FontSize',12)

ax3 = nexttile(t5);
plot(data.all.chla,data.all.ER,'.','color',rgb('lightgrey'))
hold on
plot(data.summer.chla,data.summer.ER,'.k')
ylim([-300 100])
pbaspect(ax3,[1 1 1]); % Make relative lengths of axes equal
% make a text object for the title
xl = get(gca(),'Xlim');
yl = get(gca(),'Ylim');
text(xl(1),yl(2),'Summer','VerticalAlignment','top','HorizontalAlignment','left','FontSize',12)

ax4 = nexttile(t5);
plot(data.all.chla,data.all.ER,'.','color',rgb('lightgrey'))
hold on
plot(data.fall.chla,data.fall.ER,'.k')
ylim([-300 100])
pbaspect(ax4,[1 1 1]); % Make relative lengths of axes equal
% make a text object for the title
xl = get(gca(),'Xlim');
yl = get(gca(),'Ylim');
text(xl(1),yl(2),'Fall','VerticalAlignment','top','HorizontalAlignment','left','FontSize',12)

% Control x- and y-tick labels
ax = t5.Children;
xlbl = get(ax,'XLabel');
ylbl = get(ax,'YLabel');
set([xlbl{:}],'String','');
set([ylbl{:}],'String','');
set(ax,'XTickLabel',{})
set(ax,'YTickLabel',{})
set([ax(1),ax(2)],'XTickLabelMode','auto')
set([ax(2),ax(4)],'YTickLabelMode','auto')

xlabel(t5,'Chl a (RFU)')

%----Tile 6: Turbidity-----------------------------------------------------
t6 = tiledlayout(t,2,2,'Tilespacing','none');
t6.Layout.Tile = 6;

ax1 = nexttile(t6);
plot(data.all.turbidity,data.all.ER,'.','color',rgb('lightgrey'))
hold on
plot(data.winter.turbidity,data.winter.ER,'.k')
ylim([-300 100])
pbaspect(ax1,[1 1 1]); % Make relative lengths of axes equal
% make a text object for the title
xl = get(gca(),'Xlim');
yl = get(gca(),'Ylim');
text(xl(1),yl(2),'Winter','VerticalAlignment','top','HorizontalAlignment','left','FontSize',12)

ax2 = nexttile(t6);
plot(data.all.turbidity,data.all.ER,'.','color',rgb('lightgrey'))
hold on
plot(data.spring.turbidity,data.spring.ER,'.k')
ylim([-300 100])
pbaspect(ax2,[1 1 1]); % Make relative lengths of axes equal
% make a text object for the title
xl = get(gca(),'Xlim');
yl = get(gca(),'Ylim');
text(xl(1),yl(2),'Spring','VerticalAlignment','top','HorizontalAlignment','left','FontSize',12)

ax3 = nexttile(t6);
plot(data.all.turbidity,data.all.ER,'.','color',rgb('lightgrey'))
hold on
plot(data.summer.turbidity,data.summer.ER,'.k')
ylim([-300 100])
pbaspect(ax3,[1 1 1]); % Make relative lengths of axes equal
% make a text object for the title
xl = get(gca(),'Xlim');
yl = get(gca(),'Ylim');
text(xl(1),yl(2),'Summer','VerticalAlignment','top','HorizontalAlignment','left','FontSize',12)

ax4 = nexttile(t6);
plot(data.all.turbidity,data.all.ER,'.','color',rgb('lightgrey'))
hold on
plot(data.fall.turbidity,data.fall.ER,'.k')
ylim([-300 100])
pbaspect(ax4,[1 1 1]); % Make relative lengths of axes equal
% make a text object for the title
xl = get(gca(),'Xlim');
yl = get(gca(),'Ylim');
text(xl(1),yl(2),'Fall','VerticalAlignment','top','HorizontalAlignment','left','FontSize',12)

% Control x- and y-tick labels
ax = t6.Children;
xlbl = get(ax,'XLabel');
ylbl = get(ax,'YLabel');
set([xlbl{:}],'String','');
set([ylbl{:}],'String','');
set(ax,'XTickLabel',{})
set(ax,'YTickLabel',{})
set([ax(1),ax(2)],'XTickLabelMode','auto')
set([ax(2),ax(4)],'YTickLabelMode','auto')

xlabel(t6,'Turbidity (NTU)')

ylabel(t,'ER (mmol O_2 m^{-2} d^{-1})')

%%
ax2 = nexttile;
plot(gull_dat.temperature,gull_dat.ER,'.','color',rgb('lightgrey'))
hold on
plot(gull_dat.temperature(indSpring),gull_dat.ER(indSpring),'.k')
title('Spring')

ax3 = nexttile;
plot(gull_dat.temperature,gull_dat.ER,'.','color',rgb('lightgrey'))
hold on
plot(gull_dat.temperature(indSummer),gull_dat.ER(indSummer),'.k')
title('Summer')

ax4 = nexttile;
plot(gull_dat.temperature,gull_dat.ER,'.','color',rgb('lightgrey'))
hold on
plot(gull_dat.temperature(indFall),gull_dat.ER(indFall),'.k')
title('Fall')

title(t,siteName(i),'fontsize',16)
xlabel(t,'Temperature (^oC)','fontsize',14)
ylabel(t,'ER (mmol O_2 m^{-2} d^{-1})','fontsize',14)