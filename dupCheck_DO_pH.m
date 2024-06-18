%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dupCheck_DO_pH.m
% This script compares the QC'd results from the BC and ERDC sondes for the
% recalculated DO concentration and pH and produces the "best guess" time 
% series data for each of these parameters.
%
% AUTHOR:
% Emily Chua 
% 
% DATE:
% First created: 4/26/2024
% Last updated: 5/15/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;close all;clc

rootpath = 'G:\My Drive\Postdoc\Work\SMIIL\';

prompt = {'Choose the platform'};
site = questdlg(prompt,'Platform Selection','Gull','North','South','Gull');

%====Import cleaned BC and ERDC data=======================================
cd([rootpath,'open-water-platform-data\',site,'\cleaned\movmed\'])
load([site,'-bc-cleaned'])
load([site,'-erdc-cleaned'])
dat1 = sonde1_cleaned;
dat2 = sonde2_cleaned;
% Synchronize the BC and ERDC data to a common datetime vector
dat_syn = synchronize(dat1,dat2);

clearvars dat1 dat2 sonde1_cleaned sonde2_cleaned

%====Import best-guess data for depth, S, and T============================
cd([rootpath,'open-water-platform-data\',site,'\cleaned\dupcheck\'])
load([site,'-bestGuess.mat'])

%====Import re-calculated DO concentrations================================
cd([rootpath,'open-water-platform-data\',site,'\cleaned\corrected'])
load([site,'-DOCorr.mat'])

%====Global plotting settings==============================================
dt1 = datetime('29-Jun-2021','TimeZone','UTC');     % Make all plots have same start date
dt2 = dat_syn.datetime_utc(end);
red = [0.8500 0.3250 0.0980];   % BC sonde
blue = [0 0.4470 0.7410];       % ERDC sonde
fontsize = 14;
linewidth = 1;
dotsize = 12;
circlesize = 4;

% Find indices of deployment changes
ind_dep = find(diff(bestguess.deployment.deployment) > 0);
switch site
    case 'Gull'
        dep = table([1;dat_syn.deployment_dat1(ind_dep+1)],[1;ind_dep+1]);
    case 'North'
        dep = table([2;dat_syn.deployment_dat1(ind_dep+1)],[1;ind_dep+1]);
    case 'South'
        dep = table([1;dat_syn.deployment_dat1(ind_dep+1)],[1;ind_dep+1]);
end
dep.Properties.VariableNames = {'depNum','ind'};

switch site
    case 'Gull'
        label = {'Deployment 1','Deployment 2','Deployment 5','Deployment 6',...
            'Deployment 7','Deployment 8','Deployment 9','Deployment 10',...
            'Deployment 11','Deployment 12','Deployment 13','Deployment 14',...
            'Deployment 15','Deployment 16','Deployment 17'};
    case 'North'
        label = {'Deployment 2','Deployment 6','Deployment 7','Deployment 8',...
            'Deployment 9','Deployment 10','Deployment 11','Deployment 12',...
            'Deployment 13','Deployment 14','Deployment 15','Deployment 16','Deployment 17'};
    case 'South'
        label = {'Deployment 1','Deployment 2','Deployment 4','Deployment 5',...
            'Deployment 6','Deployment 7','Deployment 8','Deployment 9',...
            'Deployment 10','Deployment 11','Deployment 12','Deployment 13',...
            'Deployment 14','Deployment 15','Deployment 16','Deployment 17'};
end

% Compare DO and pH time series
% fig = figure(1);clf
% fig.WindowState = 'maximized';
% t = tiledlayout(2,1,'TileSpacing','tight');
% 
% ax1 = nexttile;
% plot(DOcorr.bc.datetime_utc,DOcorr.bc.DOconc,'.-','Color',rgb('darkred'),'MarkerSize',dotsize,'LineWidth',linewidth,'DisplayName','BC - Re-calculated')
% hold on
% plot(DOcorr.erdc.datetime_utc,DOcorr.erdc.DOconc,'.-','Color',rgb('darkblue'),'MarkerSize',dotsize,'LineWidth',linewidth,'DisplayName','ERDC - Re-calculated')
% xline([bestguess.deployment.datetime_utc(1); bestguess.deployment.datetime_utc(ind_dep+1)],'--',label,'HandleVisibility','off')
% ylabel('DO Conc (\mumol/L)')
% legend('show')

% ax2 = nexttile;
% plot(dat_syn.datetime_utc,dat_syn.pH_dat1,'.-','Color',red,'MarkerSize',dotsize,'LineWidth',linewidth,'DisplayName','BC')
% hold on
% plot(dat_syn.datetime_utc,dat_syn.pH_dat2,'.-','Color',blue,'MarkerSize',dotsize,'LineWidth',linewidth,'DisplayName','ERDC')
% xline([bestguess.deployment.datetime_utc(1); bestguess.deployment.datetime_utc(ind_dep+1)],'--',label,'HandleVisibility','off')
% ylabel('pH')
% legend('show')
% 
% linkaxes([ax1 ax2],'x')

%====Create "best-guess" DO conc===========================================
% Find mean of DO time series 
DOcorr_syn = synchronize(DOcorr.bc,DOcorr.erdc);
DOcorr_syn.Properties.VariableNames = {'bc','erdc'};
mean_DOcorr = mean(DOcorr_syn(:,1:2),2);
DOcorr_syn.mean = mean_DOcorr.mean;
dt_utc = DOcorr_syn.datetime_utc;

% Calculate standard deviation for Deployments 11 & 12
t1 = dat_syn.datetime_utc(ind_dep(8));
t2 = DOcorr_syn.datetime_utc;
ind_d11 = interp1(t2,1:length(t2),t1,'nearest');    % Find start index for Dep 11

t1 = dat_syn.datetime_utc(ind_dep(10));
ind_d13 = interp1(t2,1:length(t2),t1,'nearest');    % Find end index for Dep 12

stdev_DO = std(DOcorr_syn(ind_d11:ind_d13,1:2),0,2,'omitmissing');
stdev.DO = mean(stdev_DO.std,'omitmissing');

% Find indices of deployment changes in synchronized, corrected DO data
t1 = dat_syn.datetime_utc(dep.ind);
t2 = DOcorr_syn.datetime_utc;
ind_depDO = interp1(t2,1:length(t2),t1,'nearest');
ind_depDO = [1;ind_depDO(2:end)];

cd([rootpath,'figures\open-water-platform\',site,'\data-qc\synchronized\DOconc'])
% Check agreement between BC and ERDC sondes
DO_bc = DOcorr_syn.bc;
DO_erdc = DOcorr_syn.erdc;
threshold = 20; % umol/L
ind_diverge = find(abs(DO_bc - DO_erdc) > threshold);

pc_exceed = height(ind_diverge)/height(DOcorr_syn)*100; % Percentage of synchronized data that exceed threshold

txt = ['Differs by >',num2str(threshold),' \mumol/L (',num2str(pc_exceed,'%.1f'),'%)'];

fig = figure(1);clf
fig.WindowState = 'maximized';
plot(dt_utc,DOcorr_syn.mean,'-','Color',rgb('magenta'),'DisplayName','Mean Value')
hold on
plot(dt_utc,DO_bc,'.','Color',rgb('darkred'),'MarkerSize',dotsize,'DisplayName','BC (recalculated)');
plot(dt_utc,DO_erdc,'.','Color',rgb('darkblue'),'MarkerSize',dotsize,'DisplayName','ERDC (recalculated)')
plot(dt_utc(ind_diverge),DO_bc(ind_diverge),'o','color',rgb('goldenrod'),'MarkerSize',circlesize,'DisplayName',txt)
xline(dt_utc(ind_depDO),'--',label,'HandleVisibility','off')
ylabel('DO concentration (\mumol/L)')
legend('show','location','best')
title([site,' - After Initial QC and Moving Median Test'])
xlim([dt1 dt2])                 % Use same x limits for comparing sites
set(gca,'FontSize',fontsize)

% Choose "better" deployment
switch site
    case 'Gull'
        d1 = DOcorr_syn.erdc(ind_depDO(1):ind_depDO(2));       % Dep 1
        d2 = DOcorr_syn.mean(ind_depDO(2)+1:ind_depDO(3));     % Dep 2
        d5 = DOcorr_syn.mean(ind_depDO(3)+1:ind_depDO(4));     % Dep 5
        d6 = DOcorr_syn.bc(ind_depDO(4)+1:ind_depDO(5));       % Dep 6
        d7 = DOcorr_syn.bc(ind_depDO(5)+1:ind_depDO(6));       % Dep 7
        d8 = DOcorr_syn.erdc(ind_depDO(6)+1:ind_depDO(7));     % Dep 8
        d9 = DOcorr_syn.mean(ind_depDO(7)+1:ind_depDO(8));     % Dep 9
        d10 = DOcorr_syn.mean(ind_depDO(8)+1:ind_depDO(9));    % Dep 10
        d11 = DOcorr_syn.mean(ind_depDO(9)+1:ind_depDO(10));   % Dep 11
        d12 = DOcorr_syn.mean(ind_depDO(10)+1:ind_depDO(11));  % Dep 12
        d13 = DOcorr_syn.erdc(ind_depDO(11)+1:ind_depDO(12));  % Dep 13
        d14 = DOcorr_syn.erdc(ind_depDO(12)+1:ind_depDO(13));  % Dep 14
        d15 = DOcorr_syn.mean(ind_depDO(13)+1:ind_depDO(14));  % Dep 15
        d16 = DOcorr_syn.mean(ind_depDO(14)+1:ind_depDO(15));  % Dep 16
        d17 = DOcorr_syn.mean(ind_depDO(15)+1:end);            % Dep 17

        DO_bestguess = [d1;d2;d5;d6;d7;d8;d9;d10;d11;d12;d13;d14;d15;d16;d17];

    case 'North'
        d2 = DOcorr_syn.erdc(ind_depDO(1):ind_depDO(2));       % Dep 2
        d6 = DOcorr_syn.bc(ind_depDO(2)+1:ind_depDO(3));       % Dep 6
        d7 = DOcorr_syn.mean(ind_depDO(3)+1:ind_depDO(4));     % Dep 7
        d8 = DOcorr_syn.mean(ind_depDO(4)+1:ind_depDO(5));     % Dep 8
        d9 = DOcorr_syn.mean(ind_depDO(5)+1:ind_depDO(6));     % Dep 9
        d10 = DOcorr_syn.mean(ind_depDO(6)+1:ind_depDO(7));    % Dep 10
        d11 = DOcorr_syn.mean(ind_depDO(7)+1:ind_depDO(8));    % Dep 11
        d12 = DOcorr_syn.mean(ind_depDO(8)+1:ind_depDO(9));    % Dep 12
        d13 = DOcorr_syn.erdc(ind_depDO(9)+1:ind_depDO(10));   % Dep 13
        d14 = DOcorr_syn.erdc(ind_depDO(10)+1:ind_depDO(11));  % Dep 14
        d15 = DOcorr_syn.mean(ind_depDO(11)+1:ind_depDO(12));  % Dep 15
        d16 = DOcorr_syn.mean(ind_depDO(12)+1:ind_depDO(13));  % Dep 16
        d17 = DOcorr_syn.mean(ind_depDO(13)+1:end);            % Dep 17

        DO_bestguess = [d2;d6;d7;d8;d9;d10;d11;d12;d13;d14;d15;d16;d17];

    case 'South'
        d1 = DOcorr_syn.mean(ind_depDO(1):ind_depDO(2));       % Dep 1
        d2 = DOcorr_syn.mean(ind_depDO(2)+1:ind_depDO(3));     % Dep 2
        d4 = DOcorr_syn.bc(ind_depDO(3)+1:ind_depDO(4));       % Dep 4
        d5 = DOcorr_syn.bc(ind_depDO(4)+1:ind_depDO(5));       % Dep 5
        % d6 = DOcorr_syn.bc(ind_depDO(5)+1:ind_depDO(6));       % Dep 6
        d6 = [DOcorr_syn.bc(ind_depDO(5)+1:29593); NaN(ind_depDO(6)-29593,1)];       % Dep 6
        d7 = DOcorr_syn.erdc(ind_depDO(6)+1:ind_depDO(7));     % Dep 7
        d8 = DOcorr_syn.mean(ind_depDO(7)+1:ind_depDO(8));     % Dep 8
        d9 = DOcorr_syn.mean(ind_depDO(8)+1:ind_depDO(9));     % Dep 9
        d10 = DOcorr_syn.bc(ind_depDO(9)+1:ind_depDO(10));     % Dep 10
        d11 = DOcorr_syn.erdc(ind_depDO(10)+1:ind_depDO(11));  % Dep 11
        d12 = DOcorr_syn.erdc(ind_depDO(11)+1:ind_depDO(12));  % Dep 12
        d13 = DOcorr_syn.erdc(ind_depDO(12)+1:ind_depDO(13));  % Dep 13
        d14 = DOcorr_syn.erdc(ind_depDO(13)+1:ind_depDO(14));  % Dep 14
        d15 = NaN(abs(ind_depDO(14) - ind_depDO(15)),1);  % Dep 15
        d16 = DOcorr_syn.mean(ind_depDO(15)+1:ind_depDO(16));  % Dep 16
        d17 = DOcorr_syn.mean(ind_depDO(16)+1:end);            % Dep 17

        DO_bestguess = [d1;d2;d4;d5;d6;d7;d8;d9;d10;d11;d12;d13;d14;d15;d16;d17];
end

clearvars d1 d2 d4 d5 d6 d7 d8 d9 d10 d11 d12 d13 d14 d15 d16 d17

DO_bestguess = table(DOcorr_syn.datetime_utc,DO_bestguess);
DO_bestguess.Properties.VariableNames = {'datetime_utc','DOconc'};
DO_bestguess = table2timetable(DO_bestguess);
DO_bestguess = rmmissing(DO_bestguess);

cd([rootpath,'figures\open-water-platform\',site,'\data-qc\synchronized\DOconc'])
fig = figure(3);clf
fig.WindowState = 'maximized';
plot(DOcorr_syn.datetime_utc,DOcorr_syn.bc,'.','Color',red,'MarkerSize',dotsize,'LineWidth',linewidth,'DisplayName','BC (recalculated)')
hold on
plot(DOcorr_syn.datetime_utc,DOcorr_syn.erdc,'.','Color',blue,'MarkerSize',dotsize,'LineWidth',linewidth,'DisplayName','ERDC (recalculated)')
plot(DOcorr_syn.datetime_utc,DOcorr_syn.mean,'.','Color',rgb('magenta'),'MarkerSize',dotsize,'LineWidth',linewidth,'DisplayName','Mean')
plot(DO_bestguess.datetime_utc,DO_bestguess.DOconc,'.k','MarkerSize',dotsize,'LineWidth',linewidth,'DisplayName','"Best Guess"')
xline(dt_utc(ind_depDO),'--',label,'HandleVisibility','off')
xlim([dt1 dt2])                 % Use same x limits for comparing sites
ylim([0 450])
ylabel('DO Conc (\mumol/L)')
title([site,' - "Best-Guess" DO Concentration'])
set(gca,'FontSize',fontsize)
legend('show','location','best')

%====Create "best-guess" pH================================================
% Find mean of pH time series 
bc = table(dat_syn.datetime_utc,dat_syn.pH_dat1);
bc.Properties.VariableNames = {'datetime_utc','pH'};
bc = table2timetable(bc);
bc = rmmissing(bc);

erdc = table(dat_syn.datetime_utc,dat_syn.pH_dat2);
erdc.Properties.VariableNames = {'datetime_utc','pH'};
erdc = table2timetable(erdc);
erdc = rmmissing(erdc);

pH_syn = synchronize(bc,erdc);
pH_syn.Properties.VariableNames = {'bc','erdc'};
mean_pH = mean(pH_syn(:,1:2),2);
pH_syn.mean = mean_pH.mean;
dt_utc = pH_syn.datetime_utc;

% Calculate standard deviation for Deployments 11 & 12
t1 = dat_syn.datetime_utc(ind_dep(8));
t2 = pH_syn.datetime_utc;
ind_d11 = interp1(t2,1:length(t2),t1,'nearest');    % Find start index for Dep 11

t1 = dat_syn.datetime_utc(ind_dep(10));
ind_d13 = interp1(t2,1:length(t2),t1,'nearest');    % Find end index for Dep 12

stdev_pH = std(pH_syn(ind_d11:ind_d13,1:2),0,2,'omitmissing');
stdev.pH = mean(stdev_pH.std,'omitmissing');

% Find indices of deployment changes in synchronized pH data
t1 = dat_syn.datetime_utc(dep.ind);
t2 = pH_syn.datetime_utc;
ind_deppH = interp1(t2,1:length(t2),t1,'nearest');
ind_deppH = [1;ind_deppH(2:end)];

cd([rootpath,'figures\open-water-platform\',site,'\data-qc\synchronized\pH'])
% Check agreement between BC and ERDC sondes
pH_bc = pH_syn.bc;
pH_erdc = pH_syn.erdc;
threshold = 0.3; 
ind_diverge = find(abs(pH_bc - pH_erdc) > threshold);

pc_exceed = height(ind_diverge)/height(pH_syn)*100; % Percentage of synchronized data that exceed threshold

txt = ['Differs by >',num2str(threshold),' (',num2str(pc_exceed,'%.1f'),'%)'];

fig = figure(3);clf
fig.WindowState = 'maximized';
plot(dt_utc,pH_syn.mean,'-','Color',rgb('magenta'),'DisplayName','Mean Value')
hold on
plot(dt_utc,pH_bc,'.','Color',red,'MarkerSize',dotsize,'DisplayName','BC');
plot(dt_utc,pH_erdc,'.','Color',blue,'MarkerSize',dotsize,'DisplayName','ERDC')
plot(dt_utc(ind_diverge),pH_bc(ind_diverge),'o','color',rgb('goldenrod'),'MarkerSize',circlesize,'DisplayName',txt)
xline(dt_utc(ind_deppH),'--',label,'HandleVisibility','off')
ylabel('pH')
legend('show','location','best')
title([site,' - After Initial QC and Moving Median Test'])
xlim([dt1 dt2])                 % Use same x limits for comparing sites
set(gca,'FontSize',fontsize)

% Choose "better" deployment
switch site
    case 'Gull'
        d1 = pH_syn.erdc(ind_deppH(1):ind_deppH(2));       % Dep 1
        d2 = pH_syn.erdc(ind_deppH(2)+1:ind_deppH(3));     % Dep 2
        d5 = pH_syn.bc(ind_deppH(3)+1:ind_deppH(4));       % Dep 5
        d6 = pH_syn.mean(ind_deppH(4)+1:ind_deppH(5));     % Dep 6
        d7 = pH_syn.mean(ind_deppH(5)+1:ind_deppH(6));     % Dep 7
        d8 = pH_syn.erdc(ind_deppH(6)+1:ind_deppH(7));     % Dep 8
        d9 = pH_syn.mean(ind_deppH(7)+1:ind_deppH(8));     % Dep 9
        d10 = pH_syn.mean(ind_deppH(8)+1:ind_deppH(9));    % Dep 10
        d11 = pH_syn.mean(ind_deppH(9)+1:ind_deppH(10));   % Dep 11
        d12 = pH_syn.mean(ind_deppH(10)+1:ind_deppH(11));  % Dep 12
        d13 = [pH_syn.erdc(ind_deppH(11)+1:69205); pH_syn.bc(69206:ind_deppH(12))];  % Dep 13
        d14 = pH_syn.erdc(ind_deppH(12)+1:ind_deppH(13));  % Dep 14
        d15 = pH_syn.mean(ind_deppH(13)+1:ind_deppH(14));  % Dep 15
        d16 = pH_syn.bc(ind_deppH(14)+1:ind_deppH(15));    % Dep 16
        d17 = pH_syn.erdc(ind_deppH(15)+1:end);            % Dep 17

        pH_bestguess = [d1;d2;d5;d6;d7;d8;d9;d10;d11;d12;d13;d14;d15;d16;d17];

    case 'North'
        d2 = pH_syn.erdc(ind_deppH(1):ind_deppH(2));       % Dep 2
        d6 = pH_syn.bc(ind_deppH(2)+1:ind_deppH(3));       % Dep 6
        d7 = pH_syn.erdc(ind_deppH(3)+1:ind_deppH(4));     % Dep 7
        d8 = pH_syn.mean(ind_deppH(4)+1:ind_deppH(5));     % Dep 8
        d9 = pH_syn.mean(ind_deppH(5)+1:ind_deppH(6));     % Dep 9
        d10 = pH_syn.mean(ind_deppH(6)+1:ind_deppH(7));    % Dep 10
        d11 = pH_syn.bc(ind_deppH(7)+1:ind_deppH(8));      % Dep 11
        d12 = pH_syn.mean(ind_deppH(8)+1:ind_deppH(9));    % Dep 12
        d13 = pH_syn.erdc(ind_deppH(9)+1:ind_deppH(10));   % Dep 13
        d14 = pH_syn.bc(ind_deppH(10)+1:ind_deppH(11));    % Dep 14
        d15 = pH_syn.mean(ind_deppH(11)+1:ind_deppH(12));  % Dep 15
        d16 = pH_syn.bc(ind_deppH(12)+1:ind_deppH(13));    % Dep 16
        d17 = pH_syn.mean(ind_deppH(13)+1:end);            % Dep 17

        pH_bestguess = [d2;d6;d7;d8;d9;d10;d11;d12;d13;d14;d15;d16;d17];
        
        % Remove "obvious" outliers
        pH_bestguess(811) = NaN;
        pH_bestguess(6422:6576) = NaN;
        
    case 'South'
        % Funky deployments
        % Dep 13, 14, and 15 -- offset adjustment
        ind1 = ind_deppH(12)-1;
        ind2 = find(~isnan(pH_syn.erdc(ind_deppH(12):end)),1);
        ind2 = ind1 + ind2;
        delta = pH_syn.erdc(ind1-1) - pH_syn.erdc(ind2);
        
        d1 = pH_syn.erdc(ind_deppH(1):ind_deppH(2));       % Dep 1
        d2 = pH_syn.erdc(ind_deppH(2)+1:ind_deppH(3));     % Dep 2
        d4 = NaN;                                          % Dep 4
        d5 = pH_syn.erdc(ind_deppH(4)+1:ind_deppH(5));     % Dep 5
        d6 = pH_syn.bc(ind_deppH(5)+1:ind_deppH(6));       % Dep 6
        d7 = pH_syn.bc(ind_deppH(6)+1:ind_deppH(7));       % Dep 7
        d8 = pH_syn.mean(ind_deppH(7)+1:ind_deppH(8));     % Dep 8
        d9 = pH_syn.mean(ind_deppH(8)+1:ind_deppH(9));     % Dep 9
        d10 = pH_syn.bc(ind_deppH(9)+1:ind_deppH(10));     % Dep 10
        d11 = pH_syn.erdc(ind_deppH(10)+1:ind_deppH(11));  % Dep 11
        d12 = pH_syn.erdc(ind_deppH(11)+1:ind_deppH(12));  % Dep 12
        d13 = pH_syn.erdc(ind_deppH(12)+1:ind_deppH(13)) + delta;  % Dep 13
        d14 = pH_syn.erdc(ind_deppH(13)+1:ind_deppH(14)) + delta;  % Dep 14
        d15 = pH_syn.erdc(ind_deppH(14)+1:ind_deppH(15)) + delta;  % Dep 15
        d16 = pH_syn.mean(ind_deppH(15)+1:ind_deppH(16));  % Dep 16
        d17 = pH_syn.mean(ind_deppH(16)+1:end);            % Dep 17

        pH_bestguess = [d1;d2;d4;d5;d6;d7;d8;d9;d10;d11;d12;d13;d14;d15;d16;d17];
end

clearvars d1 d2 d5 d6 d7 d8 d9 d10 d11 d12 d13 d14 d15 d16 d17

pH_bestguess = table(pH_syn.datetime_utc,pH_bestguess);
pH_bestguess.Properties.VariableNames = {'datetime_utc','pH'};
pH_bestguess = table2timetable(pH_bestguess);
pH_bestguess = rmmissing(pH_bestguess);

cd([rootpath,'figures\open-water-platform\',site,'\data-qc\synchronized\pH'])
fig = figure(3);clf
fig.WindowState = 'maximized';
plot(pH_syn.datetime_utc,pH_syn.bc,'.','Color',red,'MarkerSize',dotsize,'LineWidth',linewidth,'DisplayName','BC')
hold on
plot(pH_syn.datetime_utc,pH_syn.erdc,'.','Color',blue,'MarkerSize',dotsize,'LineWidth',linewidth,'DisplayName','ERDC')
plot(pH_syn.datetime_utc,pH_syn.mean,'.','Color',rgb('magenta'),'MarkerSize',dotsize,'LineWidth',linewidth,'DisplayName','Mean')
plot(pH_bestguess.datetime_utc,pH_bestguess.pH,'.k','MarkerSize',dotsize,'LineWidth',linewidth,'DisplayName','"Best Guess"')
xline(dt_utc(ind_deppH),'--',label,'HandleVisibility','off')
xlim([dt1 dt2])                 % Use same x limits for comparing sites
ylabel('pH')
title([site,' - "Best Guess" pH'])
set(gca,'FontSize',fontsize)
legend('show','location','best')

%====Save "best guess" data================================================
% Add "best-guess" DO and pH into structure
bestguess.DOconc = DO_bestguess;
bestguess.pH = pH_bestguess;

option = questdlg('Save best guess data?','Save File','Yes','No','Yes');
switch option
    case 'Yes'
        cd([rootpath,'open-water-platform-data\',site,'\cleaned\dupcheck'])
        save([site,'-bestGuess.mat'],'bestguess','stdev')
        disp('File saved!')
    case 'No'
        disp('File not saved.')
end

cd([rootpath,'figures\open-water-platform\',site,'\data-qc\synchronized\DOconc'])
