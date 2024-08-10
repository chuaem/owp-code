%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get_S_T_for_DIC_TA.m
% This script pulls out the final quality controlled salinity and temperature
% values for the open water platforms at the sample times corresponding to the
% discrete DIC and TA samples.
%
% AUTHOR:
% Emily Chua
%
% DATE:
% First created: 7/31/2024
% Last updated: 8/10/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;close all;clc

rootpath = 'G:\My Drive\Postdoc\Work\SMIIL\';

%====Import data===========================================================
cd([rootpath,'discrete-samples'])
filename = 'SMIIL DIC_TA Sample Inventory.xlsx';
inventory = readtable(filename,'Sheet','Inventory');

cd([rootpath,'open-water-platform-data\north\cleaned\final-qc'])
load('north-cleaned.mat');
tbl_all.north = finalQC;

cd([rootpath,'open-water-platform-data\gull\cleaned\final-qc'])
load('gull-cleaned.mat');
tbl_all.gull = finalQC;

cd([rootpath,'open-water-platform-data\south\cleaned\final-qc'])
load('south-cleaned.mat');
tbl_all.south = finalQC;

clearvars finalQC

% Create datetimes for all samples in local and UTC format
T  = datetime(inventory.Sampling_Time,'ConvertFrom','datenum','Format','hh:mm:ss');
datetime_local = inventory.Date_Sampled + timeofday(T);
datetime_local.TimeZone = 'America/New_York';
datetime_utc = datetime_local;
datetime_utc.TimeZone = 'UTC';

% Pull out rows corresponding to North, Gull, and South
ind.north = find(strcmp('North',inventory.Location_ID));
ind.gull = find(strcmp('Gull',inventory.Location_ID));
ind.south = find(strcmp('South',inventory.Location_ID));

fn = fieldnames(tbl_all);

for i = 1:numel(fn)
    % Find QC'd S and T values for closest matching datetimes to discrete samples for each site
    nearest_ind = interp1(tbl_all.(fn{i}).datetime_utc,1:length(tbl_all.(fn{i}).datetime_utc),datetime_utc(ind.(fn{i})),'nearest');

    % Overwrite the QC'd S and T value into the 'inventory' table
    inventory.Salinity(ind.(fn{i})) = tbl_all.(fn{i}).salinity(nearest_ind);
    inventory.Temperature(ind.(fn{i})) = tbl_all.(fn{i}).temperature(nearest_ind);
end

% Write table to Excel spreadsheet
cd([rootpath,'discrete-samples'])
writetable(inventory,filename,'Sheet','Inventory Modified')