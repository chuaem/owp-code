%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotAllTimeseries.m
% This script...
%
% AUTHOR:
% Emily Chua
%
% DATE:
% First created: 9/16/2024
% Last updated:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;close all;clc

rootpath = 'G:\My Drive\Postdoc\Work\SMIIL\';

%====Import data===========================================================
prompt = {'Choose the platform'};
site = questdlg(prompt,'Platform Selection','Gull','North','South','Gull');

% Load the sonde data
cd([rootpath,'\diel-method\owp-data\final-qc'])
load([site,'_obs.mat']);
dat_gull_10min = dat;
dat_gull_daily = retime(dat,'daily','mean');
dat_gull_monthly = retime(dat,'monthly','mean');
clear dat

% Load the MATLAB metabolism results
cd([rootpath,'diel-method\sensitivity-analysis\',site,'\wann_vary'])
load('diel_res.mat');
metab_gull = diel_dtd;
clear diel_dtd

% % Load the PAR data
% cd([rootpath,'physical-data/final-dataset'])
% load('par.mat')

%==== Combine daily sonde data and metabolic rates into one table========
daily_gull = synchronize(dat_gull_daily,metab_gull);
daily_gull = removevars(daily_gull,{'deployment','daystart_dt','dayend_dt','daylength','R_hourly','P_hourly','R_daily','P_daily','Tair'});
% daily_gull(daily_gull.datetime_utc < datetime(2021,6,29,"TimeZone","UTC"),:) = [];
% daily_gull(daily_gull.datetime_utc > datetime(2024,6,4,"TimeZone","UTC"),:) = [];

% ====Calculate seasonal means=============================================
% Metabolism data
% Get month numbers
mo = month(daily_gull.datetime_utc);
% Get year numbers
yr = year(daily_gull.datetime_utc);

% Define the seasons by month
isWinter = ismember(mo,[12 1 2]);
isSpring = ismember(mo,[3 4 5]);
isSummer = ismember(mo,[6 7 8]);
isFall = ismember(mo,[9 10 11]);

indWinter = find(isWinter);
indSpring = find(isSpring);
indSummer = find(isSummer);
indFall = find(isFall);

% Add a column for month
daily_gull.month = mo;

% Add a column for season number
daily_gull.season(indWinter) = 1;
daily_gull.season(indSpring) = 2;
daily_gull.season(indSummer) = 3;
daily_gull.season(indFall) = 4;
daily_gull.season = categorical(daily_gull.season);

% Add a column for year
% If statement so December is included in following winter
for i = 1:height(daily_gull)
    if daily_gull.month(i) == 12
        daily_gull.year(i) = categorical(yr(i) + 1);
    else
        daily_gull.year(i) = categorical(yr(i));
    end
end

% Find rows belonging to each year
ind_2021 = find(daily_gull.year == '2021');
ind_2022 = find(daily_gull.year == '2022');
ind_2023 = find(daily_gull.year == '2023');
ind_2024 = find(daily_gull.year == '2024');

% Find seasonal means
% G = groupsummary(T,groupvars,method,datavars)
seasonal_means_2021 = groupsummary(daily_gull(ind_2021,:),"season","mean",["depth","salinity","temperature","DOconc","DOsat","pH","chla","turbidity","GPP" "ER" "NEM"]);
seasonal_means_2022 = groupsummary(daily_gull(ind_2022,:),"season","mean",["depth","salinity","temperature","DOconc","DOsat","pH","chla","turbidity","GPP" "ER" "NEM"]);
seasonal_means_2023 = groupsummary(daily_gull(ind_2023,:),"season","mean",["depth","salinity","temperature","DOconc","DOsat","pH","chla","turbidity","GPP" "ER" "NEM"]);
seasonal_means_2024 = groupsummary(daily_gull(ind_2024,:),"season","mean",["depth","salinity","temperature","DOconc","DOsat","pH","chla","turbidity","GPP" "ER" "NEM"]);

seasonal_sd_2021 = groupsummary(daily_gull(ind_2021,:),"season","std",["depth","salinity","temperature","DOconc","DOsat","pH","chla","turbidity","GPP" "ER" "NEM"]);
seasonal_sd_2022 = groupsummary(daily_gull(ind_2022,:),"season","std",["depth","salinity","temperature","DOconc","DOsat","pH","chla","turbidity","GPP" "ER" "NEM"]);
seasonal_sd_2023 = groupsummary(daily_gull(ind_2023,:),"season","std",["depth","salinity","temperature","DOconc","DOsat","pH","chla","turbidity","GPP" "ER" "NEM"]);
seasonal_sd_2024 = groupsummary(daily_gull(ind_2024,:),"season","std",["depth","salinity","temperature","DOconc","DOsat","pH","chla","turbidity","GPP" "ER" "NEM"]);

seasonal_sd_2021(:,1:2) = [];
seasonal_sd_2022(:,1:2) = [];
seasonal_sd_2023(:,1:2) = [];
seasonal_sd_2024(:,1:2) = [];

seasonal_means_2021.year = repelem(categorical(2021), height(seasonal_means_2021))';
seasonal_means_2022.year = repelem(categorical(2022), height(seasonal_means_2022))';
seasonal_means_2023.year = repelem(categorical(2023), height(seasonal_means_2023))';
seasonal_means_2024.year = repelem(categorical(2024), height(seasonal_means_2024))';

dat_gull_seasonal =  [seasonal_means_2021 seasonal_sd_2021; seasonal_means_2022 seasonal_sd_2022; seasonal_means_2023 seasonal_sd_2023; seasonal_means_2024 seasonal_sd_2024];

dat_gull_seasonal = movevars(dat_gull_seasonal,"year",'After','season');

% Plot the seasonal means
xaxisLbl = {'F ''21',...
    'W ''22','Sp ''22','Su ''22','F ''22',...
    'W ''23','Sp ''23','Su ''23','F ''23',...
    'W ''24','Sp ''24','Su ''24'};

fig = figure(2);clf
fig.WindowState = 'maximized';
t1 = tiledlayout(4,2,'TileSpacing','tight');

ax1 = nexttile;
% plot(daily.north.datetime_utc,daily.north.temperature) % North color = lightseagreen
% hold on
errorbar(dat_gull_seasonal.mean_temperature,dat_gull_seasonal.std_temperature,'o-','MarkerSize',12,'LineWidth',2,'Color',rgb('goldenrod'),'DisplayName','Gull')
% plot(daily.south.datetime_utc,daily.south.temperature) % South color = mediumpurple
set(gca,'XTickLabel',[])
ylabel('T (^oC)')
legend('show','location','best')

ax2 = nexttile;
% plot(daily.north.datetime_utc,daily.north.pH)
% hold on
errorbar(dat_gull_seasonal.mean_pH,dat_gull_seasonal.std_pH,'o-','MarkerSize',12,'LineWidth',2,'Color',rgb('goldenrod'),'HandleVisibility','off')
% plot(daily.south.datetime_utc,daily.south.pH)
% ylim([7.4 8.6])
set(gca,'XTickLabel',[])
ylabel('pH')

ax3 = nexttile;
% plot(daily.north.datetime_utc,daily.north.DOsat)
% hold on
errorbar(dat_gull_seasonal.mean_DOsat,dat_gull_seasonal.std_DOsat,'o-','MarkerSize',12,'LineWidth',2,'Color',rgb('goldenrod'),'HandleVisibility','off')
% plot(daily.south.datetime_utc,daily.south.DOsat)
set(gca,'XTickLabel',[])
ylabel('DO (%sat)')

ax4 = nexttile;
% plot(daily.north.datetime_utc,daily.north.salinity)
% hold on
errorbar(dat_gull_seasonal.mean_salinity,dat_gull_seasonal.std_salinity,'o-','MarkerSize',12,'LineWidth',2,'Color',rgb('goldenrod'),'HandleVisibility','off')
% plot(daily.south.datetime_utc,daily.south.salinity)
set(gca,'XTickLabel',[])
ylabel('S (psu)')

ax5 = nexttile;
% plot(daily.north.datetime_utc,daily.north.chla,'DisplayName','North')
% hold on
errorbar(dat_gull_seasonal.mean_chla,dat_gull_seasonal.std_chla,'o-','MarkerSize',12,'LineWidth',2,'Color',rgb('goldenrod'),'HandleVisibility','off')
% plot(daily.south.datetime_utc,daily.south.chla,'DisplayName','South')
set(gca,'XTickLabel',[])
ylabel('Chl a (RFU)')

ax6 = nexttile;
% plot(daily.north.datetime_utc,daily.north.turbidity,'HandleVisibility','off')
% hold on
errorbar(dat_gull_seasonal.mean_turbidity,dat_gull_seasonal.std_turbidity,'o-','MarkerSize',12,'LineWidth',2,'Color',rgb('goldenrod'),'HandleVisibility','off')
% plot(daily.south.datetime_utc,daily.south.turbidity,'HandleVisibility','off')
set(gca,'XTickLabel',[])
ylabel('Turbidity (NTU)')

ax7 = nexttile;
% plot(daily.north.datetime_utc,daily.north.GPP,'-k')
% hold on
errorbar(dat_gull_seasonal.mean_GPP,dat_gull_seasonal.std_GPP,'o-','MarkerSize',12,'LineWidth',2,'Color',rgb('goldenrod'),'HandleVisibility','off')
% plot(daily.south.datetime_utc,daily.south.GPP,':k')
hold on
% plot(daily.north.datetime_utc,daily.north.ER,'-k')
errorbar(dat_gull_seasonal.mean_ER,dat_gull_seasonal.std_ER,'o-','MarkerSize',12,'LineWidth',2,'Color',rgb('goldenrod'),'HandleVisibility','off')
% plot(daily.south.datetime_utc,daily.south.ER,':k')
set(gca,'XTickLabel',xaxisLbl)
ylabel('ER and GPP')

ax8 = nexttile;
% plot(daily.north.datetime_utc,daily.north.NEM)
% hold on
errorbar(dat_gull_seasonal.mean_NEM,dat_gull_seasonal.std_NEM,'o-','MarkerSize',12,'LineWidth',2,'Color',rgb('goldenrod'),'HandleVisibility','off')
% plot(daily.south.datetime_utc,daily.south.NEM)
set(gca,'XTickLabel',xaxisLbl)
ylabel('NEM')

%% ====Plots of daily means==================================================
% fig = figure(1);clf
% fig.WindowState = 'maximized';
% t1 = tiledlayout(5,1,'TileSpacing','tight');
% 
% ax1 = nexttile;
% yyaxis left
% % plot(daily.north.datetime_utc,daily.north.temperature)
% % hold on
% plot(daily_gull.datetime_utc,daily_gull.temperature)
% % plot(daily.south.datetime_utc,daily.south.temperature)
% set(gca,'XTickLabel',[])
% ylabel('T (^oC)')
% 
% yyaxis right
% % plot(daily.north.datetime_utc,daily.north.pH)
% % hold on
% plot(daily_gull.datetime_utc,daily_gull.pH)
% % plot(daily.south.datetime_utc,daily.south.pH)
% ylim([7.4 8.6])
% ylabel('pH')
% 
% ax2 = nexttile;
% yyaxis left
% % plot(daily.north.datetime_utc,daily.north.DOsat)
% % hold on
% plot(daily_gull.datetime_utc,daily_gull.DOsat)
% % plot(daily.south.datetime_utc,daily.south.DOsat)
% set(gca,'XTickLabel',[])
% ylabel('DO (%sat)')
% 
% yyaxis right
% % plot(daily.north.datetime_utc,daily.north.salinity)
% % hold on
% plot(daily_gull.datetime_utc,daily_gull.salinity)
% % plot(daily.south.datetime_utc,daily.south.salinity)
% ylabel('S (psu)')
% 
% ax3 = nexttile;
% yyaxis left
% % plot(daily.north.datetime_utc,daily.north.chla,'DisplayName','North')
% % hold on
% plot(daily_gull.datetime_utc,daily_gull.chla,'DisplayName','Gull')
% % plot(daily.south.datetime_utc,daily.south.chla,'DisplayName','South')
% set(gca,'XTickLabel',[])
% ylabel('Chl a (RFU)')
% legend('show','Location','northwest')
% 
% yyaxis right
% % plot(daily.north.datetime_utc,daily.north.turbidity,'HandleVisibility','off')
% % hold on
% plot(daily_gull.datetime_utc,daily_gull.turbidity,'HandleVisibility','off')
% % plot(daily.south.datetime_utc,daily.south.turbidity,'HandleVisibility','off')
% ylabel('Turbidity (NTU)')
% 
% ax4 = nexttile;
% yyaxis left
% % plot(daily.north.datetime_utc,daily.north.GPP,'-k')
% % hold on
% plot(daily_gull.datetime_utc,daily_gull.GPP,'--')
% hold on
% % plot(daily.south.datetime_utc,daily.south.GPP,':k')
% % plot(daily.north.datetime_utc,daily.north.ER,'-k')
% plot(daily_gull.datetime_utc,daily_gull.ER,'--')
% % plot(daily.south.datetime_utc,daily.south.ER,':k')
% ylabel('ER and GPP')
% 
% yyaxis right
% % plot(daily.north.datetime_utc,daily.north.NEM)
% % hold on
% plot(daily_gull.datetime_utc,daily_gull.NEM)
% % plot(daily.south.datetime_utc,daily.south.NEM)
% ylabel('NEM')
% 
% nexttile
% yyaxis left
% plot(daily_gull.datetime_utc,daily_gull.depth)
% ylabel('Depth')
% yyaxis right
% plot(daily_gull.datetime_utc,daily_gull.par)
% ylabel('PAR')
% 
% title(t1,'Daily Means','FontSize',18,'FontWeight','bold')