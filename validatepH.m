%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% validatepH.m
% This script uses the discrete DIC/TA measurements to validate the
% AquaTroll data.
%
% AUTHOR:
% Emily Chua
%
% DATE:
% First created: 8/15/2024
% Last updated: 8/16/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;close all;clc

rootpath = 'G:\My Drive\Postdoc\Work\SMIIL\';

prompt = {'Choose the platform'};
site = questdlg(prompt,'Platform Selection','Gull','North','South','Gull');

%==========================================================================
%   Calculate pH from DIC and TA via CO2SYS
%==========================================================================

% Load the discrete DIC/TA data
cd([rootpath,'discrete-samples'])
load('DICTA_SMIIL_sample_output.mat')

% Create a table of just the OWP data
idx = find(strcmp(site,sample_output.Location_ID));

co2sys.(site) = sample_output(idx,:);

% Make into timetable
date = datetime(co2sys.(site).Date_Sampled);
timeOfDay = datetime(co2sys.(site).Sampling_Time,'ConvertFrom','datenum');
datetime_local = date + timeofday(timeOfDay);
datetime_local.TimeZone = 'America/New_York';
datetime_utc = datetime_local;
datetime_utc.TimeZone = 'UTC';

co2sys.(site).datetime_utc = datetime_utc;

co2sys.(site) = table2timetable(co2sys.(site),"RowTimes","datetime_utc");

% Calculate pH from discrete DIC/TA samples
par1type =    1; % Type "1" = "alkalinity"
par1     =    co2sys.(site).TA; % Value of the first parameter
par2type =    2; % Type "2" = "DIC"
par2     =    co2sys.(site).DIC; % Value of the second parameter
sal      =    co2sys.(site).Field_Salinity; % Salinity of the sample
tempin   =    co2sys.(site).Field_Temp; % In situ temperature
presin   =    0; % Lab pressure
tempout  =    0; % Doesn't matter here
presout  =    1; % In situ pressure
sil      =   50; % Concentration of silicate  in the sample (in umol/kg)
po4      =    2; % Concentration of phosphate in the sample (in umol/kg)
pHscale  =    1; % pH scale at which the input pH is reported ("1" means "Total Scale")
k1k2c    =    4; % Choice of H2CO3 and HCO3- dissociation constants K1 and K2 ("4" means "Mehrbach refit")
kso4c    =    1; % Choice of HSO4- dissociation constants KSO4 ("1" means "Dickson")

% Do the calculation. See CO2SYS's help for syntax and output format
A = CO2SYS(par1,par2,par1type,par2type,sal,tempin,tempout,presin,presout,sil,po4,pHscale,k1k2c,kso4c);

% Make a structure with tables containing the discrete pH values for each site
pH_discrete = table(datetime_utc,A(:,18));
pH_discrete.Properties.VariableNames = {'datetime_utc','pH'};

%==========================================================================
%   Compare discrete pH and AquaTroll pH
%==========================================================================
% Load the AquaTroll data
cd([rootpath,'open-water-platform-data\',site,'\cleaned\final-qc'])
load([site,'-cleaned.mat']);

% Global plotting settings
dt1 = datetime('29-Jun-2021','TimeZone','UTC');     % Make all plots have same start date
dt2 = finalQC.datetime_utc(end);
red = [0.8500 0.3250 0.0980];   % BC sonde
blue = [0 0.4470 0.7410];       % ERDC sonde
fontsize = 14;
linewidth = 1;
dotsize = 6;
circlesize = 6;

switch site
    case 'Gull'
        label = {'Deployment 1','Deployment 2','Deployment 5','Deployment 6',...
            'Deployment 7','Deployment 8','Deployment 9','Deployment 10',...
            'Deployment 11','Deployment 12','Deployment 13','Deployment 14',...
            'Deployment 15','Deployment 16','Deployment 17','Deployment 18'};
    case 'North'
        label = {'Deployment 2','Deployment 6','Deployment 7','Deployment 8',...
            'Deployment 9','Deployment 10','Deployment 11','Deployment 12',...
            'Deployment 13','Deployment 14','Deployment 15','Deployment 16','Deployment 17','Deployment 18'};
    case 'South'
        label = {'Deployment 1','Deployment 2','Deployment 4','Deployment 5',...
            'Deployment 6','Deployment 7','Deployment 8','Deployment 9',...
            'Deployment 10','Deployment 11','Deployment 12','Deployment 13',...
            'Deployment 14','Deployment 15','Deployment 16','Deployment 17','Deployment 18'};
end

% Find indices of deployment changes
ind_dep = find(diff(finalQC.deployment) > 0);
switch site
    case 'Gull'
        dep = table([1;finalQC.deployment(ind_dep+1)],[1;ind_dep+1]);
    case 'North'
        dep = table([2;finalQC.deployment(ind_dep+1)],[1;ind_dep+1]);
    case 'South'
        dep = table([1;finalQC.deployment(ind_dep+1)],[1;ind_dep+1]);
end
dep.Properties.VariableNames = {'depNum','ind'};

cd([rootpath,'figures\open-water-platform\',site,'\validation\dic_ta'])

% Plot AquaTroll pH data with discrete pH on top
figure(1),clf
plot(finalQC.datetime_utc,finalQC.pH,'.k','MarkerSize',8)
hold on
plot(pH_discrete.datetime_utc,pH_discrete.pH,'x');
xline(finalQC.datetime_utc(dep.ind),'--',label,'HandleVisibility','off')
xlim([dt1 dt2])    % Use same x limits
xlabel('UTC')
ylabel('pH')
title(site,'FontSize',fontsize)

% Linear regression
% Remove rows with missing pH values
finalQC_trimmed = rmmissing(finalQC,'DataVariables',"pH");

% Find closest matching indices
t1 = pH_discrete.datetime_utc;
t2 = finalQC_trimmed.datetime_utc;
ind_sonde = interp1(t2,1:length(t2),t1,'nearest','extrap');

ind_discrete = 1:1:height(pH_discrete);
% If closest matching samples are more than 25 minutes off, do not include
for i = 1:length(ind_sonde(~isnan(ind_sonde)))
    dt = between(finalQC_trimmed.datetime_utc(ind_sonde(i)),pH_discrete.datetime_utc(i));
    if abs(time(dt)) > minutes(25)
        ind_sonde(i) = NaN;
        ind_discrete(i) = NaN;
    end
end

ind_sonde = rmmissing(ind_sonde);
ind_discrete = rmmissing(ind_discrete);

% Create the model
x = finalQC_trimmed.pH(ind_sonde);
y = pH_discrete.pH(ind_discrete);
mdl = fitlm(x,y,'y~x1-1');  % Force intercept through zero; see Wilkinson notation

tbl = table(finalQC_trimmed.datetime_utc(ind_sonde),pH_discrete.datetime_utc(ind_discrete),x,y);
tbl.Properties.VariableNames = {'AquaTroll datetime','Discrete datetime','AquaTroll pH','Discrete pH'};

% Create equation string for plot
eqn = ['y = ',num2str(mdl.Coefficients.Estimate,3),'x'];
R2 = num2str(mdl.Rsquared.Ordinary,2);

fig2 = figure(2);clf
fig2.WindowState = 'maximized';
h = plot(mdl,'marker','.','markersize',20);
delete(h([3 4]))    % Delete confidence bounds on plot
hold on
plot([min(min([x,y])) max(max([x,y]))], [min(min([x,y])) max(max([x,y]))],'--k')
legend('',[eqn,newline,'R^2 = ',R2],'1:1 line','Location','southeast')
xlabel('AquaTroll pH')          
ylabel('Discrete pH')  
title(site)
daspect([1 1 1])
ylim([min(min([x,y])) max(max([x,y]))])
xlim([min(min([x,y])) max(max([x,y]))])
