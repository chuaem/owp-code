% Sample data for three time series
ts1 = [1.2, 1.3, 1.4, 1.5, 1.7, 1.8];
ts2 = [2.1, 2.3, 2.4, 2.5, 2.6, 2.7];
ts3 = [3.1, 3.3, 3.4, 3.5, 3.7, 3.8];

% Combine the data into one matrix
data = [ts1'; ts2'; ts3'];

% Create a grouping variable
group = [repmat({'TS1'}, length(ts1), 1); repmat({'TS2'}, length(ts2), 1); repmat({'TS3'}, length(ts3), 1)];

% Perform one-way ANOVA
[p, tbl, stats] = anova1(data, group);

% Display the ANOVA table
disp('ANOVA Table:');
disp(tbl);

% Perform multiple comparisons if the ANOVA is significant
if p < 0.05
    disp('Performing multiple comparisons...');
    multcompare(stats);
end