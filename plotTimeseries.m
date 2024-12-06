%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotTimeseries.m
% This script plots daily averages of the measured parameters and
% calculated metabolic rates by site and year, to compare site and 
% interannual variability. 
%
% AUTHOR:
% Emily Chua
%
% DATE:
% First created: 11/25/2024
% Last updated:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;close all;clc

rootpath = 'G:\My Drive\Postdoc\Work\SMIIL\';

%====Import final QC'd data & diel analysis results========================
fn = {'north','gull','south'};

for i = 1:length(fn)
    cd([rootpath,'open-water-platform-data\',fn{i},'\cleaned\final-qc'])
    load('finalQC.mat');
    params.(fn{i}) = finalQC;
    cd([rootpath,'diel-method\uncertainty-analysis\',fn{i}])
    load('MonteCarloResults')
    metab.(fn{i}) = diel_dtd_MC;
end
clearvars finalQC diel_dtd diel_obs

% Find daily tidal range for each site
gull_range = groupsummary(params.gull,"datetime_utc","day","range","depth");
gull_range(find(gull_range.GroupCount < 144),:) = []; % Remove days that don't have all depth measurements (full day = 144 points)
gull_range.day_datetime_utc = datetime(cellstr(gull_range.day_datetime_utc),'TimeZone','UTC');
gull_range = table2timetable(gull_range);
gull_range.GroupCount = [];
gull_range.Properties.VariableNames = {'tidal'};

north_range = groupsummary(params.north,"datetime_utc","day","range","depth");
north_range(find(north_range.GroupCount < 144),:) = []; % Remove days that don't have all depth measurements (full day = 144 points)
north_range.day_datetime_utc = datetime(cellstr(north_range.day_datetime_utc),'TimeZone','UTC');
north_range = table2timetable(north_range);
north_range.GroupCount = [];
north_range.Properties.VariableNames = {'tidal'};

south_range = groupsummary(params.south,"datetime_utc","day","range","depth");
south_range(find(south_range.GroupCount < 144),:) = []; % Remove days that don't have all depth measurements (full day = 144 points)
south_range.day_datetime_utc = datetime(cellstr(south_range.day_datetime_utc),'TimeZone','UTC');
south_range = table2timetable(south_range); 
south_range.GroupCount = [];
south_range.Properties.VariableNames = {'tidal'};

%====Import windspeed and PAR data=========================================
cd([rootpath,'physical-data\final-dataset'])
load('windspeed.mat')

% cd([rootpath,'physical-data\final-dataset'])
% load('par.mat')

%====Calculate daily means for each site===================================
gull_dailyAvg = retime(params.gull,'daily','mean');
north_dailyAvg = retime(params.north,'daily','mean');
south_dailyAvg = retime(params.south,'daily','mean');

wspd_dailyAvg = retime(era5Dat,'daily','mean');

%====Create time tables for each site======================================
%----Gull------------------------------------------------------------------
dt2 = dateshift(metab.gull.daystart_dt,'start','day');
metab.gull.daystart_dt = dt2;
daily.gull = synchronize(gull_dailyAvg,gull_range,wspd_dailyAvg,metab.gull);
daily.gull = removevars(daily.gull,{'deployment','daystart_dt','dayend_dt','daylength','R_hourly_avg','P_hourly_avg','R_daily_avg','P_daily_avg'});
daily.gull(daily.gull.datetime_utc < datetime(2021,6,29,"TimeZone","UTC"),:) = [];
daily.gull(daily.gull.datetime_utc > datetime(2024,6,4,"TimeZone","UTC"),:) = [];

% Remove timepoints with anomalous GPP and ER values
anomER = find(daily.gull.ER_avg > 0);
anomGPP = find(daily.gull.GPP_avg < 0);
daily.gull([anomER;anomGPP],:) = [];

%----North-----------------------------------------------------------------
dt2 = dateshift(metab.north.daystart_dt,'start','day');
metab.north.daystart_dt = dt2;
daily.north = synchronize(north_dailyAvg,north_range,wspd_dailyAvg,metab.north);
daily.north = removevars(daily.north,{'deployment','daystart_dt','dayend_dt','daylength','R_hourly_avg','P_hourly_avg','R_daily_avg','P_daily_avg'});
daily.north(daily.north.datetime_utc < datetime(2021,6,29,"TimeZone","UTC"),:) = [];
daily.north(daily.north.datetime_utc > datetime(2024,6,4,"TimeZone","UTC"),:) = [];

% Remove timepoints with anomalous GPP and ER values
anomER = find(daily.north.ER_avg > 0);
anomGPP = find(daily.north.GPP_avg < 0);
daily.north([anomER;anomGPP],:) = [];

%----South-----------------------------------------------------------------
dt2 = dateshift(metab.south.daystart_dt,'start','day');
metab.south.daystart_dt = dt2;
daily.south = synchronize(south_dailyAvg,south_range,wspd_dailyAvg,metab.south);
daily.south = removevars(daily.south,{'deployment','daystart_dt','dayend_dt','daylength','R_hourly_avg','P_hourly_avg','R_daily_avg','P_daily_avg'});
daily.south(daily.south.datetime_utc < datetime(2021,6,29,"TimeZone","UTC"),:) = [];
daily.south(daily.south.datetime_utc > datetime(2024,6,4,"TimeZone","UTC"),:) = [];

% Remove timepoints with anomalous GPP and ER values
anomER = find(daily.south.ER_avg > 0);
anomGPP = find(daily.south.GPP_avg < 0);
daily.south([anomER;anomGPP],:) = [];

% See triad scheme in https://files-aje-com.s3.amazonaws.com/www/row/_assets/docs/Using_Color_In_Your_Manuscript_Figures.pdf
% and https://www.rapidtables.com/convert/color/rgb-to-hex.html?r=86&g=180&b=233
% north_clr = '#56B4E9';
north_clr = '#5D3A9B';
gull_clr = '#019E73';
south_clr = '#D55E00';
% north_clr = '#017ddf';
% gull_clr = '#52d589';
% south_clr = '#b46d15';

t = datetime(2021,01,01):calyears(1):datetime(2023,01,01);
yr = year(t);

xaxisLbl = t;
xtk = 1:3;

cd([rootpath,'figures'])

%% ====Plot daily means of the parameters==================================
fig = figure(2);clf
fig.WindowState = 'maximized';
% t1 = tiledlayout(6,1,'tilespacing','tight');
t1 = tiledlayout(8,1,'TileSpacing','tight');
 
ax1 = nexttile;
ymin = -5;
ymax = 35;
plot(daily.north.datetime_utc,daily.north.temperature,'-','LineWidth',1.5,'Color',north_clr,'DisplayName','North')
hold on
plot(daily.gull.datetime_utc,daily.gull.temperature,'-.','LineWidth',1.5,'Color',gull_clr,'DisplayName','Gull')
plot(daily.south.datetime_utc,daily.south.temperature,':','LineWidth',1.5,'Color',south_clr,'DisplayName','South')
datetick('x','yyyy')
xlim([daily.gull.datetime_utc(1) daily.gull.datetime_utc(end)])
set(gca,'XTickLabel',[])
ylim([ymin ymax])
yticks([ymin (ymin+ymax)/2 ymax])
text(0.02,0.75,'Temp','Units','normalized','VerticalAlignment','bottom','FontSize',12,'FontWeight','bold')
ylabel('^oC')
set(gca,'box','off','LineWidth',1.5)
grid on
legend('show','Location','northoutside','Orientation','horizontal')

ax2 = nexttile;
ymin = 25;
ymax = 40;
plot(daily.north.datetime_utc,daily.north.salinity,'-','LineWidth',1.5,'Color',north_clr)
hold on
plot(daily.gull.datetime_utc,daily.gull.salinity,'-.','LineWidth',1.5,'Color',gull_clr)
plot(daily.south.datetime_utc,daily.south.salinity,':','LineWidth',1.5,'Color',south_clr)
datetick('x','yyyy')
xlim([daily.gull.datetime_utc(1) daily.gull.datetime_utc(end)])
set(gca,'XTickLabel',[])
ylim([ymin ymax])
yticks([ymin (ymin+ymax)/2 ymax])
text(0.02,0.75,'Sal','Units','normalized','VerticalAlignment','bottom','FontSize',12,'FontWeight','bold')
ylabel('psu')
set(gca,'box','off','LineWidth',1.5)
grid on

ax3 = nexttile;
ymin = 25;
ymax = 125;
plot(daily.north.datetime_utc,daily.north.DOsat,'-','LineWidth',1.5,'Color',north_clr)
hold on
plot(daily.gull.datetime_utc,daily.gull.DOsat,'-.','LineWidth',1.5,'Color',gull_clr)
plot(daily.south.datetime_utc,daily.south.DOsat,':','LineWidth',1.5,'Color',south_clr)
datetick('x','yyyy')
xlim([daily.gull.datetime_utc(1) daily.gull.datetime_utc(end)])
set(gca,'XTickLabel',[])
ylim([ymin ymax])
yticks([ymin (ymin+ymax)/2 ymax])
text(0.02,0.75,'DO sat','Units','normalized','VerticalAlignment','bottom','FontSize',12,'FontWeight','bold')
ylabel('%sat')
set(gca,'box','off','LineWidth',1.5)
grid on

ax4 = nexttile;
ymin = 7.3;
ymax = 8.7;
plot(daily.north.datetime_utc,daily.north.pH,'-','LineWidth',1.5,'Color',north_clr)
hold on
plot(daily.gull.datetime_utc,daily.gull.pH,'-.','LineWidth',1.5,'Color',gull_clr)
plot(daily.south.datetime_utc,daily.south.pH,':','LineWidth',1.5,'Color',south_clr)
datetick('x','yyyy')
xlim([daily.gull.datetime_utc(1) daily.gull.datetime_utc(end)])
set(gca,'XTickLabel',[])
ylim([ymin ymax])
yticks([ymin (ymin+ymax)/2 ymax])
text(0.02,0.75,'pH','Units','normalized','VerticalAlignment','bottom','FontSize',12,'FontWeight','bold')
% ylabel('-')
set(gca,'box','off','LineWidth',1.5)
grid on

ax5 = nexttile;
ymin = .5;
ymax = 2.5;
plot(daily.north.datetime_utc,daily.north.tidal,'-','LineWidth',1.5,'Color',north_clr)
hold on
plot(daily.gull.datetime_utc,daily.gull.tidal,'-.','LineWidth',1.5,'Color',gull_clr)
plot(daily.south.datetime_utc,daily.south.tidal,':','LineWidth',1.5,'Color',south_clr)
datetick('x','yyyy')
xlim([daily.gull.datetime_utc(1) daily.gull.datetime_utc(end)])
set(gca,'XTickLabel',[])
ylim([ymin ymax])
yticks([ymin (ymin+ymax)/2 ymax])
text(0.02,0.75,'Tidal','Units','normalized','VerticalAlignment','bottom','FontSize',12,'FontWeight','bold')
ylabel('m')
set(gca,'box','off','LineWidth',1.5)
grid on

ax6 = nexttile;
ymin = 0;
ymax = 15;
plot(daily.south.datetime_utc,daily.south.wspd,'-k','LineWidth',1.5)
datetick('x','yyyy')
xlim([daily.gull.datetime_utc(1) daily.gull.datetime_utc(end)])
set(gca,'XTickLabel',[])
ylim([ymin ymax])
yticks([ymin (ymin+ymax)/2 ymax])
text(0.02,0.755,'Wind speed','Units','normalized','VerticalAlignment','bottom','FontSize',12,'FontWeight','bold')
ylabel('m/s')
set(gca,'box','off','LineWidth',1.5)
grid on

ax7 = nexttile;
ymin = -750;
ymax = 450;
% ymin = -700;
% ymax = 700;
plot(daily.north.datetime_utc,daily.north.GPP_avg,'-','LineWidth',1.5,'Color',north_clr)
hold on
plot(daily.gull.datetime_utc,daily.gull.GPP_avg,'-.','LineWidth',1.5,'Color',gull_clr)
plot(daily.south.datetime_utc,daily.south.GPP_avg,':','LineWidth',1.5,'Color',south_clr)
plot(daily.north.datetime_utc,daily.north.ER_avg,'-','LineWidth',1.5,'Color',north_clr)
plot(daily.gull.datetime_utc,daily.gull.ER_avg,'--','LineWidth',1.5,'Color',gull_clr)
plot(daily.south.datetime_utc,daily.south.ER_avg,':','LineWidth',1.5,'Color',south_clr)
yline(0,'k','LineWidth',1.5)
datetick('x','yyyy')
xlim([daily.gull.datetime_utc(1) daily.gull.datetime_utc(end)])
set(gca,'XTickLabel',[])
ylim([ymin ymax])
yticks([ymin (ymin+ymax)/2 ymax])
text(0.02,0.75,'GPP','Units','normalized','VerticalAlignment','bottom','FontSize',12,'FontWeight','bold')
text(0.02,0.05,'ER','Units','normalized','VerticalAlignment','bottom','FontSize',12,'FontWeight','bold')
ylabel('mmol m^{-2} d^{-1}')
set(gca,'box','off','LineWidth',1.5)
grid on

ax8 = nexttile;
ymin = -500;
ymax = 200;
% ymax = 500;
plot(daily.north.datetime_utc,daily.north.NEM_avg,'-','LineWidth',1.5,'Color',north_clr)
hold on
plot(daily.gull.datetime_utc,daily.gull.NEM_avg,'-.','LineWidth',1.5,'Color',gull_clr)
plot(daily.south.datetime_utc,daily.south.NEM_avg,':','LineWidth',1.5,'Color',south_clr)
yline(0,'k','LineWidth',1.5)
datetick('x','yyyy')
xlim([daily.gull.datetime_utc(1) daily.gull.datetime_utc(end)])
ylim([ymin ymax])
yticks([ymin (ymin+ymax)/2 ymax])
text(0.02,0.75,'NEM','Units','normalized','VerticalAlignment','bottom','FontSize',12,'FontWeight','bold')
ylabel('mmol m^{-2} d^{-1}')
set(gca,'box','off','LineWidth',1.5)
grid on

% ax7 = nexttile;
% plot(daily.north.datetime_utc,daily.north.chla,'-','LineWidth',1.5,'Color',north_clr)
% hold on
% plot(daily.gull.datetime_utc,daily.gull.chla,'-','LineWidth',1.5,'Color',gull_clr)
% plot(daily.south.datetime_utc,daily.south.chla,'-','LineWidth',1.5,'Color',south_clr)
% datetick('x','yyyy')
% xlim([daily.gull.datetime_utc(1) daily.gull.datetime_utc(end)])
% set(gca,'XTickLabel',[])
% ylabel('Chl a (RFU)')
% 
% ax8 = nexttile;
% plot(daily.north.datetime_utc,daily.north.turbidity,'-','LineWidth',1.5,'Color',north_clr)
% hold on
% plot(daily.gull.datetime_utc,daily.gull.turbidity,'-','LineWidth',1.5,'Color',gull_clr)
% plot(daily.south.datetime_utc,daily.south.turbidity,'-','LineWidth',1.5,'Color',south_clr)
% set(gca,'XTickLabel',[])
% ylabel('Turb (NTU)')
% datetick('x','yyyy')
% xlim([daily.gull.datetime_utc(1) daily.gull.datetime_utc(end)])

set(gcf, 'Units', 'Inches', 'Position', [-9.9000   -2.6583    8.5000   10.6500], 'PaperUnits', 'Inches', 'PaperSize', [8.5, 11])

xlabel('Calendar Year')

% %% ===Plot daily means of ONLY the metabolic rates=========================
% fig = figure(2);clf
% fig.WindowState = 'maximized';
% t2 = tiledlayout(2,1,'tilespacing','tight');
% 
% ax1 = nexttile;
% plot(daily.north.datetime_utc,daily.north.GPP_avg,'-','LineWidth',1.5,'Color',north_clr)
% hold on
% plot(daily.gull.datetime_utc,daily.gull.GPP_avg,'-','LineWidth',1.5,'Color',gull_clr)
% plot(daily.south.datetime_utc,daily.south.GPP_avg,'-','LineWidth',1.5,'Color',south_clr)
% plot(daily.north.datetime_utc,daily.north.ER_avg,'-','LineWidth',1.5,'Color',north_clr)
% plot(daily.gull.datetime_utc,daily.gull.ER_avg,'-','LineWidth',1.5,'Color',gull_clr)
% plot(daily.south.datetime_utc,daily.south.ER_avg,'-','LineWidth',1.5,'Color',south_clr)
% yline(0,'k','LineWidth',1.5)
% datetick('x','yyyy')
% xlim([daily.gull.datetime_utc(1) daily.gull.datetime_utc(end)])
% set(gca,'XTickLabel',[])
% ylabel('ER & GPP')
% 
% legend('show','Location','northoutside','Orientation','horizontal')
% 
% ax2 = nexttile;
% plot(daily.north.datetime_utc,daily.north.NEM_avg,'-','LineWidth',1.5,'Color',north_clr)
% hold on
% plot(daily.gull.datetime_utc,daily.gull.NEM_avg,'-','LineWidth',1.5,'Color',gull_clr)
% plot(daily.south.datetime_utc,daily.south.NEM_avg,'-','LineWidth',1.5,'Color',south_clr)
% yline(0,'k','LineWidth',1.5)
% datetick('x','yyyy')
% xlim([daily.gull.datetime_utc(1) daily.gull.datetime_utc(end)])
% ylabel('NEM')

%% ====Bar plots of metabolic rates on daily & monthly time scales=========
%====Calculate day of year averages======================================
% Create daily averages across years for each site/sonde
dayofyear.gull = groupsummary(daily.gull,'datetime_utc','dayofyear','mean');
dayofyear.north = groupsummary(daily.north,'datetime_utc','dayofyear','mean');
dayofyear.south = groupsummary(daily.south,'datetime_utc','dayofyear','mean');

dayofyear.gull.Properties.VariableNames(1) = {'day'};
dayofyear.north.Properties.VariableNames(1) = {'day'};
dayofyear.south.Properties.VariableNames(1) = {'day'};

dayofyear_avg = table;
dayofyear_avg.day = (1:365)';
for i = 1:height(dayofyear.gull)
    dayofyear_avg.GPP(i) = mean([dayofyear.north.mean_GPP_avg(i) dayofyear.gull.mean_GPP_avg(i) dayofyear.south.mean_GPP_avg(i)],'omitmissing');
    dayofyear_avg.ER(i) = mean([dayofyear.north.mean_ER_avg(i) dayofyear.gull.mean_ER_avg(i) dayofyear.south.mean_ER_avg(i)],'omitmissing');
    dayofyear_avg.NEM(i) = mean([dayofyear.north.mean_NEM_avg(i) dayofyear.gull.mean_NEM_avg(i) dayofyear.south.mean_NEM_avg(i)],'omitmissing');
end

%====Calculate month of year averages======================================
% Create monthly bins
monthly.gull = retime(daily.gull,'monthly','mean');
monthly.north = retime(daily.north,'monthly','mean');
monthly.south = retime(daily.south,'monthly','mean');

% Create monthly averages across years for each site/sonde
monthofyear.gull = groupsummary(daily.gull,'datetime_utc','monthofyear','mean');
monthofyear.north = groupsummary(daily.north,'datetime_utc','monthofyear','mean');
monthofyear.south = groupsummary(daily.south,'datetime_utc','monthofyear','mean');

monthofyear.gull.Properties.VariableNames(1) = {'month'};
monthofyear.north.Properties.VariableNames(1) = {'month'};
monthofyear.south.Properties.VariableNames(1) = {'month'};

monthofyear_avg = table;
monthofyear_avg.month = (1:12)';
for i = 1:height(monthofyear.gull)
    monthofyear_avg.GPP(i) = mean([monthofyear.north.mean_GPP_avg(i) monthofyear.gull.mean_GPP_avg(i) monthofyear.south.mean_GPP_avg(i)]);
    monthofyear_avg.ER(i) = mean([monthofyear.north.mean_ER_avg(i) monthofyear.gull.mean_ER_avg(i) monthofyear.south.mean_ER_avg(i)]);
    monthofyear_avg.NEM(i) = mean([monthofyear.north.mean_NEM_avg(i) monthofyear.gull.mean_NEM_avg(i) monthofyear.south.mean_NEM_avg(i)]);
end

%====Create bar plots======================================================
% % Version 1
% fig = figure;clf
% t = tiledlayout(2,1);
% 
% nexttile
% bar(dayofyear_avg.day,dayofyear_avg.GPP,0.5,'k')
% hold on
% bar(dayofyear_avg.day,dayofyear_avg.ER,0.5,'k')
% bar(dayofyear_avg.day,dayofyear_avg.NEM,0.5,'FaceColor',[0.8353, 0.3686, 0])
% % ylim([-550 300])
% title('Daily')
% set(gca,'box','off','LineWidth',1.5)
% xlabel('Day of Year')
% 
% nexttile
% bar(monthofyear_avg.month,monthofyear_avg.GPP,0.5,'k')
% hold on
% bar(monthofyear_avg.month,monthofyear_avg.ER,0.5,'k')
% bar(monthofyear_avg.month,monthofyear_avg.NEM,0.5,'FaceColor',[0.8353, 0.3686, 0])
% title('Monthly')
% set(gca,'box','off','LineWidth',1.5)
% xlabel('Month of Year')
% 
% ylabel(t,'GPP, ER, and NEM (mmol O_2 m^{-2} d^{-1})','FontSize',16)
% 
% set(gcf,'Position',[-1200 -100 400 300])
% fig.Units               = 'centimeters';
% fig.Position(3)         = 18;
% fig.Position(4)         = 22;

% Version 2
% Colors from https://personal.sron.nl/~pault/
gpp_clr = '#55AA22';
er_clr = '#BB0011';

fig = figure;clf
t = tiledlayout(2,2,'TileSpacing','compact','padding','compact');

nexttile
bar(dayofyear_avg.day,dayofyear_avg.GPP,0.5,'FaceColor',gpp_clr,'EdgeColor',gpp_clr)
hold on
bar(dayofyear_avg.day,dayofyear_avg.ER,0.5,'FaceColor',er_clr,'EdgeColor',er_clr)
text(0.03,0.9,'GPP','Units','normalized','VerticalAlignment','bottom','FontSize',12,'FontWeight','bold','color',gpp_clr)
text(0.03,0.05,'ER','Units','normalized','VerticalAlignment','bottom','FontSize',12,'FontWeight','bold','color',er_clr)
set(gca,'XTickLabel',[])
ylabel('mmol O_2 m^{-2} d^{-1}','fontsize',14)
title('Daily')
set(gca,'box','off','LineWidth',1.5)

nexttile
bar(monthofyear_avg.month,monthofyear_avg.GPP,0.5,'FaceColor',gpp_clr,'EdgeColor',gpp_clr)
hold on
bar(monthofyear_avg.month,monthofyear_avg.ER,0.5,'FaceColor',er_clr,'EdgeColor',er_clr)
text(0.03,0.9,'GPP','Units','normalized','VerticalAlignment','bottom','FontSize',12,'FontWeight','bold','color',gpp_clr)
text(0.03,0.05,'ER','Units','normalized','VerticalAlignment','bottom','FontSize',12,'FontWeight','bold','color',er_clr)
set(gca,'XTickLabel',[])
title('Monthly')
set(gca,'box','off','LineWidth',1.5)

nexttile
bar(dayofyear_avg.day,dayofyear_avg.NEM,0.5,'k','EdgeColor','k')
text(0.03,0.05,'NEM','Units','normalized','VerticalAlignment','bottom','FontSize',12,'FontWeight','bold')
set(gca,'box','off','LineWidth',1.5)
ylabel('mmol O_2 m^{-2} d^{-1}','fontsize',14)
xlabel('Day of Year','fontsize',14)

nexttile
bar(monthofyear_avg.month,monthofyear_avg.NEM,0.5,'k')
text(0.03,0.05,'NEM','Units','normalized','VerticalAlignment','bottom','FontSize',12,'FontWeight','bold')
set(gca,'box','off','LineWidth',1.5)
xlabel('Month of Year','fontsize',14)

set(gcf,'Position',[-1200 -100 400 300])
fig.Units               = 'centimeters';
fig.Position(3)         = 25.2095 ;
fig.Position(4)         = 18;
