%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotDielResults.m
% This script plots the results from dielAnalyis.m
%
% AUTHOR:
% Emily Chua 
% 
% DATE:
% First created: 2/9/2024
% Last updated: 5/31/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;close all;clc

rootpath = 'G:\My Drive\Postdoc\Work\SMIIL\';

%==========================================================================
%   Import data
%==========================================================================
% Load the MATLAB detided diel analysis results for all sites
% Gull
site = 'Gull';
cd([rootpath,'\diel-method\matlab-results\final-qc\',site])
load('diel_res.mat')
diel_dtd.dayend_dt = [];
daily.gull = diel_dtd;

% North
site = 'North';
cd([rootpath,'\diel-method\matlab-results\final-qc\',site])
load('diel_res.mat')
diel_dtd.dayend_dt = [];
daily.north = diel_dtd;

% South
site = 'South';
cd([rootpath,'\diel-method\matlab-results\final-qc\',site])
load('diel_res.mat')
diel_dtd.dayend_dt = [];
daily.south = diel_dtd;

% Define colors for plots
green = "#008837";
purple = "#7b3294";
red = '#d7191c';
orange = '#fdae61';
blue = '#2c7bb6';
teal = '#5ab4ac';
brown = '#d8b365';
mauve = '#756bb1';

% Define dates for plots
start_date = datetime(2021,08,01,'TimeZone','UTC');
% end_date = datetime(2024,02,01,'TimeZone','UTC');
end_date = datetime(2024,03,05,'TimeZone','UTC');

clearvars diel_dtd diel_obs

%==========================================================================
%   Define periods to remove
%==========================================================================
% Gull
anomER = find(daily.gull.ER > 0);
anomGPP = find(daily.gull.GPP < 0);
daily.gull{[anomER;anomGPP],:} = NaN;
% ind1 = 447;
% ind_manual = ind1;
% ind_rm = [anomER; anomGPP; ind_manual];
% daily.gull{ind_rm,:} = NaN;

% North
anomER = find(daily.north.ER > 0);
anomGPP = find(daily.north.GPP < 0);
daily.north{[anomER;anomGPP],:} = NaN;
% ind1 = 248;
% ind2 = 257;
% ind_manual = [ind1;ind2];
% ind_rm = [anomER; anomGPP; ind_manual];
% daily.north{ind_rm,:} = NaN;

% South (DO data really funky for this sonde, hence all the manual removals -- see day-night_obs.fig)
anomER = find(daily.south.ER > 0);
anomGPP = find(daily.south.GPP < 0);
daily.south{[anomER;anomGPP],:} = NaN;
% ind1 = (109:111)';
% ind2 = 261;
% ind3 = 274;
% ind4 = (309:319)';
% ind5 = 362;
% ind6 = (452:454)';
% ind7 = (598:657)';
% ind8 = (682:720)';
% ind9 = 782;
% ind_manual = [ind1;ind2;ind3;ind4;ind5;ind6;ind7;ind8;ind9];
% ind_rm = [anomER; anomGPP; ind_manual];
% daily.south{ind_rm,:} = NaN;

%==========================================================================
%   Do calculations
%==========================================================================
% Create monthly bins
monthly.gull = retime(daily.gull,'monthly','mean');
monthly.north = retime(daily.north,'monthly','mean');
monthly.south = retime(daily.south,'monthly','mean');

% Create monthly averages across years for each site/sonde
monthlyavg.gull = groupsummary(daily.gull,'daystart_dt','monthofyear','mean');
monthlyavg.north = groupsummary(daily.north,'daystart_dt','monthofyear','mean');
monthlyavg.south = groupsummary(daily.south,'daystart_dt','monthofyear','mean');

monthlyavg.gull.Properties.VariableNames(1) = {'month'};
monthlyavg.north.Properties.VariableNames(1) = {'month'};
monthlyavg.south.Properties.VariableNames(1) = {'month'};

monthlystd.gull = groupsummary(daily.gull,'daystart_dt','monthofyear','std');
monthlystd.north = groupsummary(daily.north,'daystart_dt','monthofyear','std');
monthlystd.south = groupsummary(daily.south,'daystart_dt','monthofyear','std');

monthlystd.gull.Properties.VariableNames(1) = {'month'};
monthlystd.north.Properties.VariableNames(1) = {'month'};
monthlystd.south.Properties.VariableNames(1) = {'month'};

%==========================================================================
%  Calculate overall rates (average of all values for each site)
%==========================================================================
% Units: [mmol O2 m-2 d-1]
overallNEM_gull = mean(daily.gull.NEM,'omitmissing');
overallNEM_north = mean(daily.north.NEM,'omitmissing');
overallNEM_south = mean(daily.south.NEM,'omitmissing');

% Convert to [mol O2 m-2 y-1]
overallNEM_gull = overallNEM_gull*365/1000;
overallNEM_north = overallNEM_north*365/1000;
overallNEM_south  = overallNEM_south*365/1000;

%==========================================================================
%   Daily and monthly bar plots of detided data
%==========================================================================
cd('G:\My Drive\Postdoc\Work\SMIIL\figures\diel-analysis\final-qc')

%====Gull==================================================================
site = 'Gull';
fig1 = figure(1);clf
t = tiledlayout(2,1,'TileSpacing','compact');
fig1.WindowState = 'maximized';

% Daily plots
ax1 = nexttile;
yyaxis left
    b1 = bar(daily.gull.daystart_dt,daily.gull.GPP);
    hold on
    b2 = bar(daily.gull.daystart_dt,daily.gull.ER);
    b1(1).FaceColor = green;
    b2(1).FaceColor = purple;
    ylim([-300 300])
yyaxis right
    plot(daily.gull.daystart_dt,daily.gull.NEM,'.-','Color','k','LineWidth',1,'MarkerSize',6)
    ylim([-100 100])
% Customize the x-axis format
xlim([start_date end_date])
xticks = get(ax1,'XTick');
xticksMid = mean([xticks(1:end-1); xticks(2:end)],1);
dt1 = datetime(2021,10,01,'TimeZone','UTC');
xticksNew = sort([dt1 xticks xticksMid]);
set(ax1,'XTick',xticksNew,'FontSize',15)
xtickformat('MMM yyyy')
% Make the y-axes black
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
set(gca,'FontSize',15)
% Add legend
legend({'GPP','ER','NEM'},'FontSize',15,'Location',[0.8527,0.7806,0.06028,0.1028])

% Monthly plot
ax2 = nexttile;
yyaxis left
    b1 = bar(monthly.gull.daystart_dt,monthly.gull.GPP);
    hold on
    b2 = bar(monthly.gull.daystart_dt,monthly.gull.ER);
    b1(1).FaceColor = green;
    b2(1).FaceColor = purple;
    ylim([-300 300])
yyaxis right
    plot(monthly.gull.daystart_dt,monthly.gull.NEM,'.-','Color','k','LineWidth',2,'MarkerSize',6)
    ylim([-100 100])
% Customize the x-axis format
xlim([start_date end_date])
set(ax2,'XTick',xticksNew)
xtickformat('MMM yyyy')
% Make the y-axes black
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
set(gca,'FontSize',15)

linkaxes([ax1 ax2],'x')

% Make one common y-label for each side and customize
ax = axes(fig1);
han = gca;
han.Visible = 'off';

% Left y-label
yyaxis(ax, 'left');
ylabel('GPP and ER (mmol O_2 m^{-2} d^{-1})','FontSize',14,'Color','k')
han.YLabel.Visible = 'on';
% Move the left y-label farther from the plots
yl = get(gca,'ylabel');
pyl = get(yl,'position');
pyl(1) = 1.2*pyl(1);
set(yl,'position',pyl,'fontsize',20)

% Right y-label
yyaxis(ax, 'right');
ylabel('NEM (mmol O_2 m^{-2} d^{-1})','FontSize',14,'Color','k')
han.YLabel.Visible = 'on';
% Move the right y-label farther from the plots
yr = get(gca,'ylabel');
pyr = get(yr,'position');
pyr(1) = .97*pyr(1);
set(yr,'position',pyr,'FontSize',20)

title(t,site,'FontSize',36,'FontName','Helvetica','FontWeight','bold')

%====North==================================================================
site = 'North';
fig2 = figure(2);clf
t = tiledlayout(2,1,'TileSpacing','compact');
fig2.WindowState = 'maximized';

% Daily plot
ax1 = nexttile;
yyaxis left
    b1 = bar(daily.north.daystart_dt,daily.north.GPP);
    hold on
    b2 = bar(daily.north.daystart_dt,daily.north.ER);
    b1(1).FaceColor = green;
    b2(1).FaceColor = purple;
    ylim([-300 300])
yyaxis right
    plot(daily.north.daystart_dt,daily.north.NEM,'.-','Color','k','LineWidth',1,'MarkerSize',6)
    ylim([-100 100])
% Customize the x-axis format
xlim([start_date end_date])
xticks = get(ax1,'XTick');
xticksMid = mean([xticks(1:end-1); xticks(2:end)],1);
dt1 = datetime(2021,10,01,'TimeZone','UTC');
xticksNew = sort([dt1 xticks xticksMid]);
set(ax1,'XTick',xticksNew,'FontSize',15)
xtickformat('MMM yyyy')
% Make the y-axes black
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
% Add legend
legend({'GPP','ER','NEM'},'FontSize',14,'Location',[0.8527,0.7806,0.06028,0.1028])

% Monthly plot
ax2 = nexttile;
yyaxis left
    b1 = bar(monthly.north.daystart_dt,monthly.north.GPP);
    hold on
    b2 = bar(monthly.north.daystart_dt,monthly.north.ER);
    b1(1).FaceColor = green;
    b2(1).FaceColor = purple;
    ylim([-300 300])
yyaxis right
    plot(monthly.north.daystart_dt,monthly.north.NEM,'.-','Color','k','LineWidth',2,'MarkerSize',6)
    ylim([-100 100])
% Customize the x-axis format
xlim([start_date end_date])
set(ax2,'XTick',xticksNew)
xtickformat('MMM yyyy')
% Make the y-axes black
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
set(gca,'FontSize',15)

linkaxes([ax1 ax2],'x')

% Make one common y-label for each side and customize
ax = axes(fig2);
han = gca;
han.Visible = 'off';

% Left y-label
yyaxis(ax, 'left');
ylabel('GPP and ER (mmol O_2 m^{-2} d^{-1})','FontSize',14,'Color','k')
han.YLabel.Visible = 'on';
% Move the left y-label farther from the plots
yl = get(gca,'ylabel');
pyl = get(yl,'position');
pyl(1) = 1.2*pyl(1);
set(yl,'position',pyl,'FontSize',20)

% Right y-label
yyaxis(ax, 'right');
ylabel('NEM (mmol O_2 m^{-2} d^{-1})','FontSize',14,'Color','k')
han.YLabel.Visible = 'on';
% Move the right y-label farther from the plots
yr = get(gca,'ylabel');
pyr = get(yr,'position');
pyr(1) = .97*pyr(1);
set(yr,'position',pyr,'FontSize',20)

title(t,site,'FontSize',36,'FontName','Helvetica','FontWeight','bold')

%====South==================================================================
site = 'South';
fig3 = figure(3);clf
t = tiledlayout(2,1,'TileSpacing','compact');
fig3.WindowState = 'maximized';

% Daily plot
ax1 = nexttile;
yyaxis left
    b1 = bar(daily.south.daystart_dt,daily.south.GPP);
    hold on
    b2 = bar(daily.south.daystart_dt,daily.south.ER);
    b1(1).FaceColor = green;
    b2(1).FaceColor = purple;
    ylim([-300 300])
yyaxis right
    plot(daily.south.daystart_dt,daily.south.NEM,'.-','Color','k','LineWidth',1,'MarkerSize',6)
    ylim([-100 100])
% Customize the x-axis format
xlim([start_date end_date])
xticks = get(ax1,'XTick');
xticksMid = mean([xticks(1:end-1); xticks(2:end)],1);
dt1 = datetime(2021,10,01,'TimeZone','UTC');
xticksNew = sort([dt1 xticks xticksMid]);
set(ax1,'XTick',xticksNew,'FontSize',15)
xtickformat('MMM yyyy')
% Make the y-axes black
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
% Add legend
legend({'GPP','ER','NEM'},'FontSize',14,'Location',[0.8527,0.7806,0.06028,0.1028])

% Monthly plot
ax2 = nexttile;
yyaxis left
    b1 = bar(monthly.south.daystart_dt,monthly.south.GPP);
    hold on
    b2 = bar(monthly.south.daystart_dt,monthly.south.ER);
    b1(1).FaceColor = green;
    b2(1).FaceColor = purple;
    ylim([-300 300])
yyaxis right
    plot(monthly.south.daystart_dt,monthly.south.NEM,'.-','Color','k','LineWidth',2,'MarkerSize',6)
    ylim([-100 100])
% Customize the x-axis format
xlim([start_date end_date])
set(ax2,'XTick',xticksNew)
xtickformat('MMM yyyy')
% Make the y-axes black
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
set(gca,'FontSize',15)

linkaxes([ax1 ax2],'x')

% Make one common y-label for each side and customize
ax = axes(fig3);
han = gca;
han.Visible = 'off';

% Left y-label
yyaxis(ax, 'left');
ylabel('GPP and ER (mmol O_2 m^{-2} d^{-1})','FontSize',14,'Color','k')
han.YLabel.Visible = 'on';
% Move the left y-label farther from the plots
yl = get(gca,'ylabel');
pyl = get(yl,'position');
pyl(1) = 1.2*pyl(1);
set(yl,'position',pyl,'fontsize',20)

% Right y-label
yyaxis(ax, 'right');
ylabel('NEM (mmol O_2 m^{-2} d^{-1})','FontSize',14,'Color','k')
han.YLabel.Visible = 'on';
% Move the right y-label farther from the plots
yr = get(gca,'ylabel');
pyr = get(yr,'position');
pyr(1) = .97*pyr(1);
set(yr,'position',pyr,'FontSize',20)

title(t,site,'FontSize',36,'FontName','Helvetica','FontWeight','bold')

%==========================================================================
%   Compare seasonality in NEM across sites -- without errorbars
%==========================================================================
NEMavg = [monthlyavg.north.mean_NEM'; monthlyavg.gull.mean_NEM'; monthlyavg.south.mean_NEM']';
month = monthlyavg.gull.month;

fig4 = figure(4);clf
fig4.WindowState = 'maximized';
b = bar(month,NEMavg,1);
b(1).FaceColor = teal;
b(2).FaceColor = brown;
b(3).FaceColor = mauve;
set(gca,'FontSize',15)
legend('North','Gull','South','location','southeast','Location',[0.8415,0.8062,0.0639,0.0964])
ylabel('NEM (mmol O_2 m^{-2} d^{-1})','FontSize',20,'Color','k')
xlabel('Month','Fontsize',20)
title('All Sites','FontSize',36,'FontName','Helvetica','FontWeight','bold')


% % Compare seasonality in NEM across sites -- with errorbars
% % Errorbars -- % https://www.mathworks.com/matlabcentral/answers/1582304-putting-error-bars-on-bar-plot?s_tid=srchtitle
% NEMavg = [monthlyavg.north.mean_NEM'; monthlyavg.gull.mean_NEM'; monthlyavg.south.mean_NEM']';
% NEMstd = [monthlystd.north.std_NEM'; monthlystd.gull.std_NEM'; monthlystd.south.std_NEM']';
% month = monthlyavg.gull.month;
% 
% fig5 = figure(5);clf
% fig5.WindowState = 'maximized';
% b = bar(month,NEMavg);
% b(1).FaceColor = teal;
% b(2).FaceColor = brown;
% b(3).FaceColor = mauve;
% hold on
% 
% for k = 1:numel(b)
%     xtips = b(k).XEndPoints;
%     ytips = b(k).YEndPoints;
%     errorbar(xtips,ytips,NEMstd(:,k),'k','LineStyle','none','LineWidth',1);
% end
% % Replot the bars so they're on top of the errorbars
% c = bar(month,NEMavg,1);
% hold off
% 
% c(1).FaceColor = teal;
% c(2).FaceColor = brown;
% c(3).FaceColor = mauve;
% 
% legend('North','Gull','South','location','southeast')
% ylabel('NEM (mmol O_2 m^{-2} d^{-1})','FontSize',20,'Color','k')
% xlabel('Month','Fontsize',20)

%==========================================================================
%   Option to save figures
%==========================================================================
option = questdlg('Save plots as .png and .fig?','Save plots','Yes','No','Yes');

cd([rootpath,'figures\diel-analysis\final-qc\'])
switch option
    case 'Yes'
        cd([rootpath,'figures\diel-analysis\final-qc\gull\matlab-results'])
        saveas(fig1,'daily_monthly_dtd.fig')
        saveas(fig1,'daily_monthly_dtd.png')
        
        cd([rootpath,'figures\diel-analysis\final-qc\north\matlab-results'])
        saveas(fig2,'daily_monthly_dtd.fig')
        saveas(fig2,'daily_monthly_dtd.png')
        
        cd([rootpath,'figures\diel-analysis\final-qc\south\matlab-results'])
        saveas(fig3,'daily_monthly_dtd.fig')
        saveas(fig3,'daily_monthly_dtd.png')

        cd([rootpath,'figures\diel-analysis\final-qc\site-comparison'])
        saveas(fig4,'monthlyavg_allsites.fig')
        saveas(fig4,'monthlyavg_allsites.png')
        
        % cd('site-comparison')
        % saveas(fig5,'monthlyavg_allsites_err.fig')
        % saveas(fig5,'monthlyavg_allsites_err.png')

        disp('Plots saved!')
    
    case 'No'
        disp('Plots not saved.')
end