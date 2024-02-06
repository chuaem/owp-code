%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dataQC_gull.m
% This script performs data quality control on the adjusted merged sonde
% data for Gull.
%
% Code that requires manual input is commented with "INPUTS".
%
% AUTHOR:
% Emily Chua 
% 
% DATE:
% First created: 11/9/2023
% Last updated: 2/3/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;close all;clc

site = 'gull';

rootpath = 'G:\My Drive\Postdoc\Work\SMIIL\';

cd([rootpath,'open-water-platform-data\',site,'\adjusted\merged'])

load(['alldeps-',site,'-adj.mat'])

prompt = {'Choose which sonde to QC'};
answer = questdlg(prompt,'Sonde Selection','Sonde 1','Sonde 2','Cancel','Cancel');

switch answer
    case 'Sonde 1'
        dat = sonde1_all;
    case 'Sonde 2'
        dat = sonde2_all;
end

clear sonde1_all sonde2_all

%====Global plotting settings==============================================
dt1 = datetime('29-Jun-2021','TimeZone','UTC');     % Make all plots have same start date
dt2 = dat.datetime_utc(end);
red = [0.8500 0.3250 0.0980];   % BC sonde
blue = [0 0.4470 0.7410];       % ERDC sonde
fontsize = 14;
linewidth = 1;
dotsize = 2;
circlesize = 6;

label = {'Deployment 1','Deployment 2','Deployment 5','Deployment 6',...
    'Deployment 7','Deployment 8','Deployment 9','Deployment 10',...
    'Deployment 11','Deployment 12','Deployment 13','Deployment 14',...
    'Deployment 15','Deployment 16'};

% Find indices of deployment changes
ind_dep = find(diff(dat.deployment) > 0);

%====Create a table of flags for each sonde================================
% Start with everything as passing (flag = 1)
flags = ones(height(dat),width(dat));
flags = array2table(flags);
switch answer
    case 'Sonde 1'
        flags.Properties.VariableNames = {'deployment' 'datetime_utc' 'datetime_local' ...
            'actual_cond' 'specific_cond' 'salinity' 'resistivity' 'density' ...
            'temperature' 'barometric_p' 'p' 'depth' 'TDS' 'DO_conc' 'DO_sat' 'pO2' ...
            'pH' 'pH_raw' 'ORP' 'chla' 'nitrate' 'external_voltage' 'battery_capacity'};
        sondename = 'BC';
        cd([rootpath,'\figures\open-water-platform-figures\gull\data-qc\bc'])
    case 'Sonde 2'
        flags.Properties.VariableNames = {'deployment' 'datetime_utc' 'datetime_local' ...
            'actual_cond' 'specific_cond' 'salinity' 'resistivity' 'density' ...
            'temperature' 'barometric_p' 'p' 'depth' 'TDS' 'DO_conc' 'DO_sat' 'pO2' ...
            'pH' 'pH_raw' 'ORP' 'turbidity' 'external_voltage' 'battery_capacity'};
        sondename = 'ERDC';
        cd([rootpath,'\figures\open-water-platform-figures\gull\data-qc\erdc'])
end

%====(0) Manual removals (see notes in Collab Lab Notebook)================
% Plot original depth data to define points to remove manually
% fig = figure(1);clf
% fig.WindowState = 'maximized';
% h0 = plot(dat.datetime_utc,dat.depth,'.k','markersize',6);
% xline([dat.datetime_utc(1); dat.datetime_utc(ind_dep+1)],'--',label);
% ylabel('Depth (m)')
% title(['QC Tests: Gull ',sondename,' Sonde - Depth'])
% xlim([dt1 dt2])                 % Use same x limits for comparing sites
% set(gca,'FontSize',fontsize)
% xline([dat.datetime_utc(1); dat.datetime_utc(ind_dep+1)],'--',label);

% INPUTS
switch answer
    case 'Sonde 1'
        % Out-of-water points
        oow1 = (34474:34476)';      % Dep 5 --> 6
        oow2 = (51466:51486)';      % Dep 6 --> 7
        oow3 = (60821:60829)';      % Dep 7 --> 8
        oow4 = (67605:67643)';      % Dep 8 --> 9
        oow5 = (71503:71528)';      % Dep 9 --> 10
        oow6 = (78899);             % Dep 10 --> 11
        oow7 = (87962:87965)';      % Dep 11 --> 12
        oow8 = (104299:104326)';    % Dep 12 --> 13
        oow9 = (113499:113501)';    % Dep 13 --> 14
        oow10 = (122403:122405)';   % Dep 14 --> 15
        oow11 = (130652:130661)';   % Dep 15 --> 16
        ind_oow = [oow1;oow2;oow3;oow4;oow5;oow6;oow7;oow8;oow9;oow10;oow11];

        % Isolated points during large gap from 7/30/21 (during Dep 1) to start of Dep 2
        ind_manual1 = (7420:12590)';
        ind_manual2 = find(dat.depth(1:ind_dep(2)) < 0);
        ind_manual_global = unique([ind_manual1;ind_manual2]);
    
    case 'Sonde 2'
    % Out-of-water points
        oow1 = (34473:34478)'; % Dep 5 --> 6
        oow2 = (51063:51081)'; % Dep 6 --> 7
        oow3 = (60425:60429)'; % Dep 7 --> 8
        oow4 = (67205:67233)'; % Dep 8 --> 9
        oow5 = (71093:71100)'; % Dep 9 --> 10
        oow6 = (78470:78477)'; % Dep 10 --> 11
        oow7 = (87534:87540)'; % Dep 11 --> 12
        oow8 = (103874:103903)'; % Dep 12 --> 13
        oow9 = (113075:113077)'; % Dep 13 --> 14
        oow10 = (121978:121982)'; % Dep 14 --> 15
        oow11 = (130230:130231)'; % Dep 15 --> 16
        ind_oow = [oow1;oow2;oow3;oow4;oow5;oow6;oow7;oow8;oow9;oow10;oow11];

        % Isolated points during large gap from 7/30/21 (during Dep 1) to start of Dep 2
        ind_manual1 = (7420:12589)';
        ind_manual2 = find(dat.depth(1:ind_dep(2)) < 0);
        ind_manual_global = unique([ind_manual1;ind_manual2]);
end
clearvars oow1 oow2 oow3 oow4 oow5 oow6 oow7 oow8 oow9 oow10 oow11

% DEPTH
depth_orig = dat.depth;  % Preserve original depth data for plotting

%====(1) Gross Range Test==================================================
% INPUTS
d_low = 0;      % Lower limit (m) - flag all negative values as FAIL
d_high = 5;     % Upper limit (m) - not sure what should be

ind_low = find(dat.depth < d_low);
ind_high = find(dat.depth > d_high);
% Don't duplicate points already defined to be manually removed
ind_low = setdiff(ind_low,ind_manual_global);   
ind_high = setdiff(ind_high,ind_manual_global);
ind_grossRange = [ind_low;ind_high];
ind_grossRange = setdiff(ind_grossRange,ind_oow);

ind_manual = [ind_manual_global;ind_oow];

%====(2) Spike Test========================================================
% INPUTS
spike_threshold = 0.1;    % FAIL threshold (m)

d_depth = zeros(length(dat.depth),1);

for i = 2:length(dat.depth)-1
    % Only check points if both they and their neighbors have not already been flagged for manual removal
    if ~ismember(i,ind_manual) && ~ismember(i+1,ind_manual) && ~ismember(i-1,ind_manual)
        ref = (dat.depth(i+1)+dat.depth(i-1))/2; % Calculate the average of i+1 and i-1
        d_depth(i) = dat.depth(i) - ref; % Calculate the difference between i and ref
    end
end
% If i - ref is greater than high threshold, mark as a spike
ind_spike = find(abs(d_depth) >= spike_threshold);
[C,ig,is] = intersect(ind_grossRange,ind_spike); % Don't duplicate points already removed in Gross Range Test
ind_spike(is) = [];

% Plot histogram of depth differences
fig = figure(2);clf
fig.WindowState = 'maximized';
histogram(d_depth)
title(['Spike Test Histogram: Gull ',sondename,' Sonde - Depth'])
xlabel('d_{depth} (m)')
set(gca,'YScale','log','FontSize',fontsize)

% Plot original data with flagged points
fig = figure(3);clf
fig.WindowState = 'maximized';
h0 = plot(dat.datetime_utc,dat.depth,'.k','markersize',dotsize);
hold on
h1 = plot(dat.datetime_utc(ind_manual_global),dat.depth(ind_manual_global),'og','markersize',circlesize);
h2 = plot(dat.datetime_utc(ind_oow),dat.depth(ind_oow),'oc','markersize',circlesize);
h3 = plot(dat.datetime_utc(ind_grossRange),dat.depth(ind_grossRange),'or','markersize',circlesize);
h4 = plot(dat.datetime_utc(ind_spike),dat.depth(ind_spike),'om','markersize',circlesize);
xline([dat.datetime_utc(1); dat.datetime_utc(ind_dep+1)],'--',label);
hold off
ylabel('Depth (m)')
legend([h0 h1 h2 h3 h4],'Original Data','Manual Removal','Out of Water','Gross Range Test','Spike Test','location','best')
title(['QC Tests: Gull ',sondename,' Sonde - Depth'])
xlim([dt1 dt2])                 % Use same x limits for comparing sites
set(gca,'FontSize',fontsize)

%====Clean data============================================================
% Discard points that failed Gross Range and Spike Tests and do Manual Removals
dat.depth(ind_manual_global) = NaN;
dat.depth(ind_grossRange) = NaN;
dat.depth(ind_spike) = NaN;
dat.depth(ind_oow) = NaN;

% Flag all NaNs as FAIL
ind_fail = find(isnan(dat.depth));
flags.depth(ind_fail) = 4;

% Create structure to save flagged depth indices
depths_flagged = struct('ind_manual',ind_manual_global,'ind_grossRange',ind_grossRange,'ind_spike',ind_spike,'ind_oow',ind_oow,'d_high',d_high,'d_low',d_low,'spike_threshold',spike_threshold);

clearvars ind_high ind_low ind_grossRange ind_spike C ig is

%===Clean pressure data based on cleaned depth data========================
dat.p = mean(dat.density,'omitmissing').*1000*9.81.*dat.depth/6894.76;
flags.p(ind_fail) = 4;

% TEMPERATURE
T_orig = dat.temperature;  % Preserve original T data for plotting

% Plot the temperature data with the depth flags
% fig=figure;clf
% fig.WindowState = 'maximized';
% h0 = plot(dat.datetime_utc,T_orig,'.k');
% hold on
% h1 = plot(dat.datetime_utc(depths_flagged_all),T_orig(depths_flagged_all),'or');
% xline([dat.datetime_utc(1); dat.datetime_utc(ind_dep+1)],'--',label);
% hold off
% legend([h0 h1],'Original Data','Depth Flags','location','best')
% ylabel('Temperature (^oC)')
% title('QC Tests: Gull BC Sonde - Temperature')
% xlim([dt1 dt2])                 % Use same x limits for comparing sites 
% set(gca,'FontSize',FontSize)

%====(1) Gross Range Test==================================================
% INPUTS -- could make fancier by changing limits based on season??
T_low = -2;      % Lower limit (oC) - Freezing pt of SW at 35 ppt (http://www.csgnetwork.com/h2ofreezecalc.html)
T_high = 35;     % Upper limit (oC)

ind_low = find(dat.temperature < T_low);
ind_high = find(dat.temperature > T_high);
ind_low = setdiff(ind_low,ind_manual_global);
ind_high = setdiff(ind_high,ind_manual_global);
ind_grossRange = [ind_low;ind_high];

%====(2) Spike Test========================================================
% INPUTS
spike_threshold = 1;    % FAIL threshold (deg C)

d_T = zeros(length(dat.temperature),1);

for i = 2:length(dat.temperature)-1
    % Only check points if both they and their neighbors have not already been flagged for manual removal
    if ~ismember(i,ind_manual) && ~ismember(i+1,ind_manual) && ~ismember(i-1,ind_manual)
        ref = (dat.temperature(i+1)+dat.temperature(i-1))/2; % Calculate the average of i+1 and i-1
        d_T(i) = dat.temperature(i) - ref; % Calculate the difference between i and ref
    end
end
% If i - ref is greater than high threshold, mark as a spike
ind_spike = find(abs(d_T) >= spike_threshold);
[C,ig,is] = intersect(ind_grossRange,ind_spike); % Don't duplicate points already removed in Gross Range Test
ind_spike(is) = [];

% Plot histogram of temperature differences
fig = figure;clf
fig.WindowState = 'maximized';
histogram(d_T)
title(['Spike Test Histogram: Gull ',sondename,' Sonde - Temperature'])
xlabel('d_{temperature} (^oC)')
set(gca,'YScale','log','FontSize',fontsize)

% Plot original data with flagged points
fig = figure;clf
fig.WindowState = 'maximized';
h0 = plot(dat.datetime_utc,dat.temperature,'.k','MarkerSize',dotsize);
hold on
h1 = plot(dat.datetime_utc(ind_manual_global),dat.temperature(ind_manual_global),'og','MarkerSize',circlesize);
h2 = plot(dat.datetime_utc(ind_grossRange),dat.temperature(ind_grossRange),'or','MarkerSize',circlesize);
h3 = plot(dat.datetime_utc(ind_spike),dat.temperature(ind_spike),'om','MarkerSize',circlesize);
h4 = plot(dat.datetime_utc(ind_oow),dat.temperature(ind_oow),'oc','MarkerSize',circlesize);
xline([dat.datetime_utc(1); dat.datetime_utc(ind_dep+1)],'--',label);
hold off
ylabel('Temperature (^oC)')
legend([h0 h1 h2 h3 h4],'Original Data','Manual Removal','Gross Range Test','Spike Test','Out of Water','location','best')
title(['QC Tests: Gull ',sondename,' Sonde - Temperature'])
xlim([dt1 dt2])                 % Use same x limits for comparing sites 
set(gca,'FontSize',fontsize)

%====Clean data============================================================
% Discard points that failed Gross Range and Spike Tests and do Manual Removals
dat.temperature(ind_grossRange) = NaN;
dat.temperature(ind_manual_global) = NaN;
dat.temperature(ind_spike) = NaN;
dat.temperature(ind_oow) = NaN;

% Flag all NaNs as FAIL
ind_fail = find(isnan(dat.temperature));
flags.temperature(ind_fail) = 4;

% Create structure to save flagged temperature indices
T_flagged = struct('ind_manual',ind_manual_global,'ind_grossRange',ind_grossRange,'ind_spike',ind_spike,'ind_oow',ind_oow,'d_high',d_high,'d_low',d_low,'spike_threshold',spike_threshold);

clearvars ind_high ind_low ind_grossRange ind_spike C ig is

% DO CONCENTRATION
DOconc_orig = dat.DO_conc;  % Preserve original DO concentration data for plotting

%====(1) Gross Range Test==================================================
% INPUTS
DO_low = 0;      % Lower limit (umol/L)
DO_high = 500;     % Upper limit (umol/L); based on Winkler data and O2 solubility given S & T

ind_low = find(dat.DO_conc < DO_low);
ind_high = find(dat.DO_conc > DO_high);
ind_low = setdiff(ind_low,ind_manual_global);
ind_high = setdiff(ind_high,ind_manual_global);
ind_grossRange = [ind_low;ind_high];

%====(2) Spike Test========================================================
% INPUTS
spike_threshold = 50;    % FAIL threshold (umol/L)

d_DO = zeros(length(dat.DO_conc),1);

for i = 2:length(dat.DO_conc)-1
    % Only check points if both they and their neighbors have not already been flagged for manual removal
    if ~ismember(i,ind_manual) && ~ismember(i+1,ind_manual) && ~ismember(i-1,ind_manual)
        ref = (dat.DO_conc(i+1)+dat.DO_conc(i-1))/2; % Calculate the average of i+1 and i-1
        d_DO(i) = dat.DO_conc(i) - ref; % Calculate the difference between i and ref
    end
end
% If i - ref is greater than high threshold, mark as a spike
ind_spike = find(abs(d_DO) >= spike_threshold);
[C,ig,is] = intersect(ind_grossRange,ind_spike); % Don't duplicate points already removed in Gross Range Test
ind_spike(is) = [];

% Plot histogram of DO concentration differences
fig = figure;clf
fig.WindowState = 'maximized';
histogram(d_DO)
title(['Spike Test Histogram: Gull ',sondename,' DO Concentration'])
xlabel('d_{DO conc} (\mumol/L)')
set(gca,'YScale','log','FontSize',fontsize)

% Plot original data with flagged points
fig = figure;clf
fig.WindowState = 'maximized';
h0 = plot(dat.datetime_utc,dat.DO_conc,'.k','MarkerSize',dotsize);
hold on
h1 = plot(dat.datetime_utc(ind_manual_global),dat.DO_conc(ind_manual_global),'og','MarkerSize',circlesize);
h2 = plot(dat.datetime_utc(ind_grossRange),dat.DO_conc(ind_grossRange),'or','MarkerSize',circlesize);
h3 = plot(dat.datetime_utc(ind_spike),dat.DO_conc(ind_spike),'om','MarkerSize',circlesize);
h4 = plot(dat.datetime_utc(ind_oow),dat.DO_conc(ind_oow),'oc','MarkerSize',circlesize);
xline([dat.datetime_utc(1); dat.datetime_utc(ind_dep+1)],'--',label);
hold off
ylabel('DO Concentration (\mumol/L)')
legend([h0 h1 h2 h3 h4],'Original Data','Manual Removal','Gross Range Test','Spike Test','Out of Water','location','best')
title(['QC Tests: Gull ',sondename,' Sonde - DO Concentration'])
xlim([dt1 dt2])                 % Use same x limits for comparing sites--09876543 
set(gca,'FontSize',fontsize)

%====Clean data============================================================
% Discard points that failed Gross Range and Spike Tests and do Manual Removals
dat.DO_conc(ind_manual_global) = NaN;
dat.DO_conc(ind_grossRange) = NaN;
dat.DO_conc(ind_spike) = NaN;
dat.DO_conc(ind_oow) = NaN;

% Flag all NaNs as FAIL
ind_fail = find(isnan(dat.DO_conc));
flags.DO_conc(ind_fail) = 4;

% Create structure to save flagged DO concentration indices
DOconc_flagged = struct('ind_manual',ind_manual_global,'ind_grossRange',ind_grossRange,'ind_spike',ind_spike,'ind_oow',ind_oow,'d_high',d_high,'d_low',d_low,'spike_threshold',spike_threshold);

clearvars ind_high ind_low ind_grossRange ind_spike C ig is

% SALINITY
S_orig = dat.salinity;  % Preserve original S data for plotting

%====(1) Gross Range Test==================================================
% INPUTS
S_low = 0;      % Lower limit (psu) - flag all negative values as FAIL
S_high = 50;     % Upper limit (psu) - not sure what should be

ind_low = find(dat.salinity < S_low);
ind_high = find(dat.salinity > S_high);
% Don't duplicate points already defined to be manually removed
ind_low = setdiff(ind_low,ind_manual_global);   
ind_high = setdiff(ind_high,ind_manual_global);
ind_grossRange = [ind_low;ind_high];
ind_grossRange = setdiff(ind_grossRange,ind_oow);

ind_manual = [ind_manual_global;ind_oow];

%====(2) Spike Test========================================================
% INPUTS
spike_threshold = 10;    % FAIL threshold (psu)

d_salinity = zeros(length(dat.salinity),1);

for i = 2:length(dat.salinity)-1
    % Only check points if both they and their neighbors have not already been flagged for manual removal
    if ~ismember(i,ind_manual) && ~ismember(i+1,ind_manual) && ~ismember(i-1,ind_manual)
        ref = (dat.salinity(i+1)+dat.salinity(i-1))/2; % Calculate the average of i+1 and i-1
        d_salinity(i) = dat.salinity(i) - ref; % Calculate the difference between i and ref
    end
end
% If i - ref is greater than high threshold, mark as a spike
ind_spike = find(abs(d_salinity) >= spike_threshold);
[C,ig,is] = intersect(ind_grossRange,ind_spike); % Don't duplicate points already removed in Gross Range Test
ind_spike(is) = [];

% Plot histogram of salinity differences
fig = figure;clf
fig.WindowState = 'maximized';
histogram(d_salinity)
title(['Spike Test Histogram: Gull ',sondename,' Salinity'])
xlabel('d_{salinity} (psu)')
set(gca,'YScale','log','FontSize',fontsize)

% Plot original data with flagged points
fig = figure;clf
fig.WindowState = 'maximized';
h0 = plot(dat.datetime_utc,dat.salinity,'.k','MarkerSize',dotsize);
hold on
h1 = plot(dat.datetime_utc(ind_manual_global),dat.salinity(ind_manual_global),'og','MarkerSize',circlesize);
h2 = plot(dat.datetime_utc(ind_grossRange),dat.salinity(ind_grossRange),'or','MarkerSize',circlesize);
h3 = plot(dat.datetime_utc(ind_spike),dat.salinity(ind_spike),'om','MarkerSize',circlesize);
h4 = plot(dat.datetime_utc(ind_oow),dat.salinity(ind_oow),'oc','MarkerSize',circlesize);
xline([dat.datetime_utc(1); dat.datetime_utc(ind_dep+1)],'--',label)
hold off
ylabel('Salinity (psu)')
legend([h0 h1 h2 h3 h4],'Original Data','Manual Removal','Gross Range Test','Spike Test','Out of Water','location','best')
title(['QC Tests: Gull ',sondename,' Sonde - Salinity'])
xlim([dt1 dt2])                 % Use same x limits for comparing sites
set(gca,'FontSize',fontsize)

%====Clean data============================================================
% Discard points that failed Gross Range and Spike Tests and do Manual Removals
dat.salinity(ind_manual_global) = NaN;
dat.salinity(ind_grossRange) = NaN;
dat.salinity(ind_spike) = NaN;
dat.salinity(ind_oow) = NaN;

% Flag all NaNs as FAIL
ind_fail = find(isnan(dat.salinity));
flags.salinity(ind_fail) = 4;

% Create structure to save flagged salinity indices
S_flagged = struct('ind_manual',ind_manual_global,'ind_grossRange',ind_grossRange,'ind_spike',ind_spike,'ind_oow',ind_oow,'S_high',S_high,'S_low',S_low,'spike_threshold',spike_threshold);

clearvars ind_high ind_low ind_grossRange ind_spike C ig is

%====Convert table to timetable============================================
dat_TT = table2timetable(dat,'RowTimes','datetime_utc');
newTimeStep = minutes(10);
dat_TT = retime(dat_TT,'regular','mean','TimeStep',newTimeStep); % Calculate mean of values in each time bin

dat_TT.datetime_utc = dat_TT.datetime_utc + newTimeStep/2;
dat_TT.datetime_local = dat_TT.datetime_utc;
dat_TT.datetime_local.TimeZone = 'America/New_York';

%====Plot the cleaned data=================================================
cd([rootpath,'\figures\open-water-platform-figures\gull\data-qc\',sondename])

% Depth
fig = figure;clf
fig.WindowState = 'maximized';
% h0 = plot(dat.datetime_utc,depth_orig,'.','Color',red);
% hold on
h1 = plot(dat_TT.datetime_utc,dat_TT.depth,'.k','MarkerSize',dotsize);
xline([dat.datetime_utc(1); dat.datetime_utc(ind_dep+1)],'--',label);
% hold off
% legend([h0 h1],'Original','Retimed','location','best')
ylabel('Depth (m)')
title(['Gull ',sondename,' Sonde - Cleaned & Retimed Depth Data'])
xlim([dt1 dt2])                 % Use same x limits for comparing sites
ylim([-.5 3])
set(gca,'FontSize',fontsize)

% Temperature
fig = figure;clf
fig.WindowState = 'maximized';
% h0 = plot(dat.datetime_utc,T_orig,'.','Color',red);
% hold on
h1 = plot(dat_TT.datetime_utc,dat_TT.temperature,'.k','MarkerSize',dotsize);
xline([dat.datetime_utc(1); dat.datetime_utc(ind_dep+1)],'--',label);
% hold off
% legend([h0 h1],'Original','Retimed','location','best')
ylabel('Temperature (^oC)')
title(['Gull ',sondename,' Sonde - Cleaned & Retimed Temperature Data'])
xlim([dt1 dt2])                 % Use same x limits for comparing sites
set(gca,'FontSize',fontsize)

% DO
fig = figure;clf
fig.WindowState = 'maximized';
% h0 = plot(dat.datetime_utc,DOconc_orig,'.','Color',red);
% hold on
h1 = plot(dat_TT.datetime_utc,dat_TT.DO_conc,'.k','MarkerSize',dotsize);
xline([dat.datetime_utc(1); dat.datetime_utc(ind_dep+1)],'--',label);
% hold off
% legend([h0 h1],'Original','Retimed','location','best')
ylabel('DO Concentration (\mumol/L)')
title(['Gull ',sondename,' Sonde - Cleaned & Retimed DO Concentration Data'])
xlim([dt1 dt2])                 % Use same x limits for comparing sites
set(gca,'FontSize',fontsize)

% Salinity
fig = figure;clf
fig.WindowState = 'maximized';
% h0 = plot(dat.datetime_utc,S_orig,'.','Color',red);
% hold on
h1 = plot(dat_TT.datetime_utc,dat_TT.salinity,'.k','MarkerSize',dotsize);
xline([dat.datetime_utc(1); dat.datetime_utc(ind_dep+1)],'--',label);
% hold off
% legend([h0 h1],'Original','Retimed','location','best')
ylabel('Salinity (psu)')
title(['Gull ',sondename,' Sonde - Cleaned & Retimed Salinity Data'])
xlim([dt1 dt2])                 % Use same x limits for comparing sites
ylim([-.5 50])
set(gca,'FontSize',fontsize)

%====Save the cleaned data=================================================
option = questdlg('Save cleaned data?','Save File','Yes','No','Yes');
switch option
    case 'Yes'
        cd([rootpath,'open-water-platform-data\gull\cleaned'])
        switch answer
            case 'Sonde 1'
                sonde1_cleaned = dat_TT;
                save('gull-bc-cleaned.mat','sonde1_cleaned','flags')
            case 'Sonde 2'
                sonde2_cleaned = dat_TT;
                save('gull-erdc-cleaned.mat','sonde2_cleaned','flags')
        end
        disp('File saved!')
    case 'No'
        disp('File not saved.')
end
        