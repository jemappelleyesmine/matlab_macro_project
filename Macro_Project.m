%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% MACROECONOMICS PROJECT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Link to our github repository: https://github.com/jemappelleyesmine/matlab_macro_project

cd '/Users/yesminehachana/Desktop/Classes/Dauphine/2nd Semester/Macroeconomics MATLAB Project/Repository'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% QUARTERLY DATA ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 1. Load Quarterly GDP Datasets

% Load US and France data (full period)
US_GDP = readtable('GDP_US_FRED.xlsx');
France_GDP = readtable('GDP_France_ECB.xlsx');

% Load China data (2007 onward)
China_GDP = readtable('GDP_China_NBS.xlsx');

% Display first few rows to verify
disp('US GDP Data:'); disp(US_GDP(1:5,:));
disp('France GDP Data:'); disp(France_GDP(1:5,:));
disp('China GDP Data:'); disp(China_GDP(1:5,:));

%% 2. Data Cleaning & Alignment

% Rename columns for consistency
US_GDP.Properties.VariableNames = {'Quarter', 'US_GDP'};
France_GDP.Properties.VariableNames = {'Quarter', 'France_GDP'};
China_GDP.Properties.VariableNames = {'Quarter', 'China_GDP'};

% Convert "Quarter" column to datetime format (YYYY-MM-DD)
US_GDP.Quarter = datetime(US_GDP.Quarter, 'InputFormat', 'yyyy-MM-dd');
France_GDP.Quarter = datetime(France_GDP.Quarter, 'InputFormat', 'yyyy-MM-dd');
China_GDP.Quarter = datetime(China_GDP.Quarter, 'InputFormat', 'yyyy-MM-dd');

% Find the starting and ending dates
start_date_US_France = min(US_GDP.Quarter); % Earliest available date
start_date_China = datetime(2007,1,1); % China's data starts in 2007
end_date = max(US_GDP.Quarter); % Latest available date

% Filter US & France (Keep Full Period)
US_France_GDP = innerjoin(US_GDP, France_GDP, 'Keys', 'Quarter');

% Filter US, France & China (Start from 2007)
US_France_China_GDP = innerjoin(US_France_GDP(US_France_GDP.Quarter >= start_date_China, :), China_GDP, 'Keys', 'Quarter');

% Display merged datasets
disp('US-France Full GDP Dataset:');
disp(US_France_GDP(1:5,:));

disp('US-France-China (2007 Onward) GDP Dataset:');
disp(US_France_China_GDP(1:5,:));

%% 3. Convert GDP to Logarithms

% Convert US-France (Full Period)
log_gdp_US_France = log(table2array(US_France_GDP(:,2:3)));

% Convert US-France-China (2007 Onward)
log_gdp_US_France_China = log(table2array(US_France_China_GDP(:,2:4)));

% Store log-transformed GDP in structured tables
log_gdp_table_US_France = table(US_France_GDP.Quarter, log_gdp_US_France(:,1), log_gdp_US_France(:,2), ...
    'VariableNames', {'Quarter', 'Log_US_GDP', 'Log_France_GDP'});

log_gdp_table_US_France_China = table(US_France_China_GDP.Quarter, log_gdp_US_France_China(:,1), log_gdp_US_France_China(:,2), log_gdp_US_France_China(:,3), ...
    'VariableNames', {'Quarter', 'Log_US_GDP', 'Log_France_GDP', 'Log_China_GDP'});

% Display first few rows to verify
disp(log_gdp_table_US_France(1:5,:));
disp(log_gdp_table_US_France_China(1:5,:));

%% 4. Apply HP Filter (λ = 1600 for Quarterly Data)

lambda = 1600; % Standard HP filter value for quarterly data

% Apply the HP filter for US-France (Full Period)
[cycle_USA, trend_USA] = hpfilter(log_gdp_US_France(:,1), lambda);
[cycle_FRA, trend_FRA] = hpfilter(log_gdp_US_France(:,2), lambda);

% Apply the HP filter for US-France-China (2007 Onward)
[cycle_USA_China, trend_USA_China] = hpfilter(log_gdp_US_France_China(:,1), lambda);
[cycle_FRA_China, trend_FRA_China] = hpfilter(log_gdp_US_France_China(:,2), lambda);
[cycle_CHN, trend_CHN] = hpfilter(log_gdp_US_France_China(:,3), lambda);

% Store results in tables
business_cycle_table_US_France = table(US_France_GDP.Quarter, cycle_USA, cycle_FRA, ...
    'VariableNames', {'Quarter', 'US_Cycle', 'France_Cycle'});

business_cycle_table_US_France_China = table(US_France_China_GDP.Quarter, cycle_USA_China, cycle_FRA_China, cycle_CHN, ...
    'VariableNames', {'Quarter', 'US_Cycle', 'France_Cycle', 'China_Cycle'});

% Display first few rows
disp(business_cycle_table_US_France(1:5,:));
disp(business_cycle_table_US_France_China(1:5,:));

%% 5. Plot Business Cycle - US & France (Full Period)

figure;
plot(US_France_GDP.Quarter, cycle_USA, 'b', 'LineWidth', 1.5); hold on;
plot(US_France_GDP.Quarter, cycle_FRA, 'g', 'LineWidth', 1.5);
xlabel('Year');
ylabel('Business Cycle Component');
title('Quarterly GDP Business Cycle - US & France');
legend('USA', 'France');
grid on;

%Graph 1: Business Cycle of USA & France (1950 - Present):

%The USA (blue) and France (green) cycles show a clear co-movement, 
% especially in the post-1980 period.

%Some divergence occurs in earlier decades, but their fluctuations 
% remain broadly synchronized.

%The 2008 financial crisis and COVID-19 recession are visible as sharp 
% downturns in both countries.

%% 6. Plot Business Cycle - US, France & China (2007 Onward)

figure;
plot(US_France_China_GDP.Quarter, cycle_USA_China, 'b', 'LineWidth', 1.5); hold on;
plot(US_France_China_GDP.Quarter, cycle_FRA_China, 'g', 'LineWidth', 1.5);
plot(US_France_China_GDP.Quarter, cycle_CHN, 'r', 'LineWidth', 1.5);
xlabel('Year');
ylabel('Business Cycle Component');
title('Quarterly GDP Business Cycle - US, France & China (2007 Onward)');
legend('USA', 'France', 'China');
grid on;

%Graph 2: Business Cycle of USA, France, and China (2007 - Present):
%The USA and France continue to move together, while China (red) 
% exhibits a much more volatile cycle.

%China’s fluctuations appear less synchronized with the USA and France.

%The 2020 COVID-19 recession affected all three economies significantly, 
% but China’s recovery seems more volatile.

%% 7. Compute Business Cycle Statistics

% Extract business cycle components
us_cycle = business_cycle_table_US_France.US_Cycle;
france_cycle = business_cycle_table_US_France.France_Cycle;

% Extract China's cycle separately (2007 onward)
china_cycle = business_cycle_table_US_France_China.China_Cycle;

% Compute Volatility (Standard Deviation)
volatility_us = std(us_cycle);
volatility_france = std(france_cycle);
volatility_china = std(china_cycle);

% Compute Persistence (First-Order Autocorrelation)
persistence_us = autocorr(us_cycle, 1);
persistence_france = autocorr(france_cycle, 1);
persistence_china = autocorr(china_cycle, 1);

% Compute Co-movement (Correlation between Countries)
corr_us_france = corrcoef(us_cycle, france_cycle);

% China only has data from 2007 onward, match the period
corr_china_us = corrcoef(china_cycle, business_cycle_table_US_France_China.US_Cycle);
corr_china_france = corrcoef(china_cycle, business_cycle_table_US_France_China.France_Cycle);

% Store results in a table
disp('Business Cycle Statistics:')
business_cycle_stats = table(["USA"; "France"; "China"], ...
    [volatility_us; volatility_france; volatility_china], ...
    [persistence_us(2); persistence_france(2); persistence_china(2)], ...
    [corr_us_france(1,2); corr_china_france(1,2); corr_china_us(1,2)], ...
    'VariableNames', {'Country', 'Volatility', 'Persistence', 'Correlation'});

disp(business_cycle_stats);

%Volatility (Standard Deviation):

%China (0.07759) has the highest volatility, indicating that its GDP 
% fluctuations are more extreme compared to the USA and France.

%The USA (0.01735) and France (0.01590) have relatively lower volatility,
% suggesting that their business cycles exhibit more stability.

%The higher volatility in China could be explained by its rapid economic 
% transitions, structural shifts, and exposure to external shocks.

%Persistence (First-Order Autocorrelation):

%The USA (0.78014) has the most persistent business cycle fluctuations. 
% This means that its economic expansions and contractions tend to last 
% longer.

%France (0.57052) shows moderate persistence, suggesting that its business 
% cycle deviations do not last as long as the USA’s.

%China (0.08627) has very low persistence, meaning that its economic 
% fluctuations are short-lived and less predictable.

%Correlation (Co-Movement):

%USA and France (0.46985) have the strongest correlation, implying that 
% their economies move together relatively closely. This is likely due to 
% their high level of trade and financial integration.

%France and China (0.25548) have a weaker correlation, suggesting that 
% France’s business cycle is somewhat aligned with China’s but not 
% strongly.

%China and the USA (0.17751) show the weakest correlation, indicating that 
% their business cycles are largely independent, likely due to different 
% economic structures, trade policies, and domestic demand dynamics.

%% 8. Performing Lagged and Lead Correlations

max_lag = 12; % Define max lag to check (quarters)

% Initialize matrices to store results
corr_lags = zeros(2*max_lag+1, 3); % Rows for -12 to +12 lags

for lag = -max_lag:max_lag
    if lag < 0
        % Leading: Shift second country forward
        corr_lags(lag+max_lag+1, 1) = corr(china_cycle(1:end+lag), business_cycle_table_US_France_China.US_Cycle(-lag+1:end)); 
        corr_lags(lag+max_lag+1, 2) = corr(china_cycle(1:end+lag), business_cycle_table_US_France_China.France_Cycle(-lag+1:end)); 
        corr_lags(lag+max_lag+1, 3) = corr(us_cycle(1:end+lag), france_cycle(-lag+1:end)); 
    elseif lag > 0
        % Lagging: Shift first country forward
        corr_lags(lag+max_lag+1, 1) = corr(china_cycle(lag+1:end), business_cycle_table_US_France_China.US_Cycle(1:end-lag)); 
        corr_lags(lag+max_lag+1, 2) = corr(china_cycle(lag+1:end), business_cycle_table_US_France_China.France_Cycle(1:end-lag)); 
        corr_lags(lag+max_lag+1, 3) = corr(us_cycle(lag+1:end), france_cycle(1:end-lag)); 
    else
        % No shift (original correlation)
        corr_lags(lag+max_lag+1, 1) = corr(china_cycle, business_cycle_table_US_France_China.US_Cycle);
        corr_lags(lag+max_lag+1, 2) = corr(china_cycle, business_cycle_table_US_France_China.France_Cycle);
        corr_lags(lag+max_lag+1, 3) = corr(us_cycle, france_cycle);
    end
end

% Store results in a table
lag_labels = (-max_lag:max_lag)';
lead_lag_table = table(lag_labels, corr_lags(:,1), corr_lags(:,2), corr_lags(:,3), ...
    'VariableNames', {'Lag', 'China_US_Corr', 'China_France_Corr', 'US_France_Corr'});

disp('Lead-Lag Correlations:')
disp(lead_lag_table)

% Plot Lead-Lag Correlations
figure;
plot(lag_labels, corr_lags(:,1), '-or', 'LineWidth', 1.5); hold on;
plot(lag_labels, corr_lags(:,2), '-ob', 'LineWidth', 1.5);
plot(lag_labels, corr_lags(:,3), '-og', 'LineWidth', 1.5);
xlabel('Lag (quarters)');
ylabel('Correlation');
title('Lead-Lag Correlation of Business Cycles');
legend('China-USA', 'China-France', 'USA-France');
grid on;

%China and the USA
%The correlation is generally weak across different lags.

%The highest correlation (0.35881) appears at one quarter lag, suggesting 
% that China’s business cycle leads the US cycle slightly in some cases.

%Negative correlations at longer leads and lags suggest that their economic
% cycles do not move in sync over time.

%China and France
%The correlation is strongest at one quarter lag (0.39845), meaning China’s
% cycle slightly precedes France’s.

%The correlation is weak at other time points, showing that their business 
% cycles are not strongly aligned in general.

%USA and France
%The highest correlation (0.46985) is at lag 0, meaning their cycles move 
% in sync at the same time.

%The correlation remains positive for a few quarters before weakening, 
% reinforcing the idea that the USA and France have interdependent 
% economies.

%Graph 3: Lead-Lag Correlation:
%The USA and France have the strongest correlation at lag 0, meaning their
% economies fluctuate together without delay.

%China’s business cycle is weakly correlated with the USA and France, 
% with its impact appearing more strongly after a short delay 
% (1-2 quarters).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% DESEASONALIZED QUARTERLY DATA ANALYSIS %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%We used the moving average method to remove seasonal fluctuations 
% from the log-GDP data. Specifically, we applied a centered moving 
% average over a 4-quarter window to capture and smooth out seasonal 
% patterns. Moving averages are easy to implement which is why we tried
% this method.

%The moving average computes the average of adjacent data points over 
% a defined window (in our case, 4 quarters).

%This helps smooth out short-term seasonal variations while retaining 
% long-term trends and cycles.

%We subtracted this smoothed seasonal component from the original data 
% to obtain the deseasonalized GDP.

%% 1. Remove Seasonality using Moving Average Method

% Define moving average window (4 quarters for seasonal adjustment)
window = 4;

% Compute moving average (centered)
mov_avg_US = movmean(log_gdp_US_France(:,1), window, 'Endpoints', 'discard');
mov_avg_France = movmean(log_gdp_US_France(:,2), window, 'Endpoints', 'discard');

mov_avg_US_China = movmean(log_gdp_US_France_China(:,1), window, 'Endpoints', 'discard');
mov_avg_France_China = movmean(log_gdp_US_France_China(:,2), window, 'Endpoints', 'discard');
mov_avg_China = movmean(log_gdp_US_France_China(:,3), window, 'Endpoints', 'discard');

% Compute seasonal component
seasonal_US = log_gdp_US_France(1:length(mov_avg_US),1) - mov_avg_US;
seasonal_France = log_gdp_US_France(1:length(mov_avg_France),2) - mov_avg_France;

seasonal_US_China = log_gdp_US_France_China(1:length(mov_avg_US_China),1) - mov_avg_US_China;
seasonal_France_China = log_gdp_US_France_China(1:length(mov_avg_France_China),2) - mov_avg_France_China;
seasonal_China = log_gdp_US_France_China(1:length(mov_avg_China),3) - mov_avg_China;

% Compute seasonally adjusted log GDP
log_gdp_adj_US = log_gdp_US_France(1:length(mov_avg_US),1) - seasonal_US;
log_gdp_adj_France = log_gdp_US_France(1:length(mov_avg_France),2) - seasonal_France;

log_gdp_adj_US_China = log_gdp_US_France_China(1:length(mov_avg_US_China),1) - seasonal_US_China;
log_gdp_adj_France_China = log_gdp_US_France_China(1:length(mov_avg_France_China),2) - seasonal_France_China;
log_gdp_adj_China = log_gdp_US_France_China(1:length(mov_avg_China),3) - seasonal_China;

%% 2. Apply HP Filter to Seasonally Adjusted Data

% Apply the HP filter for US-France (Full Period)
[cycle_USA_adj, trend_USA_adj] = hpfilter(log_gdp_adj_US, lambda);
[cycle_FRA_adj, trend_FRA_adj] = hpfilter(log_gdp_adj_France, lambda);

% Apply the HP filter for US-France-China (2007 Onward)
[cycle_USA_China_adj, trend_USA_China_adj] = hpfilter(log_gdp_adj_US_China, lambda);
[cycle_FRA_China_adj, trend_FRA_China_adj] = hpfilter(log_gdp_adj_France_China, lambda);
[cycle_CHN_adj, trend_CHN_adj] = hpfilter(log_gdp_adj_China, lambda);

% Store results in new tables
business_cycle_table_US_France_adj = table(US_France_GDP.Quarter(1:length(cycle_USA_adj)), ...
    cycle_USA_adj, cycle_FRA_adj, ...
    'VariableNames', {'Quarter', 'US_Cycle', 'France_Cycle'});

business_cycle_table_US_France_China_adj = table(US_France_China_GDP.Quarter(1:length(cycle_USA_China_adj)), ...
    cycle_USA_China_adj, cycle_FRA_China_adj, cycle_CHN_adj, ...
    'VariableNames', {'Quarter', 'US_Cycle', 'France_Cycle', 'China_Cycle'});

%% 3. Compute Business Cycle Statistics for Deseasonalized Data

% Extract deseasonalized business cycle components
us_cycle_adj = cycle_USA_adj;
france_cycle_adj = cycle_FRA_adj;
china_cycle_adj = cycle_CHN_adj;

% Compute Volatility (Standard Deviation)
volatility_us_adj = std(us_cycle_adj);
volatility_france_adj = std(france_cycle_adj);
volatility_china_adj = std(china_cycle_adj);

% Compute Persistence (First-Order Autocorrelation)
persistence_us_adj = autocorr(us_cycle_adj, 1);
persistence_france_adj = autocorr(france_cycle_adj, 1);
persistence_china_adj = autocorr(china_cycle_adj, 1);

% Compute Co-movement (Correlation between Countries)
corr_us_france_adj = corrcoef(us_cycle_adj, france_cycle_adj);
corr_china_us_adj = corrcoef(china_cycle_adj, cycle_USA_China_adj);
corr_china_france_adj = corrcoef(china_cycle_adj, cycle_FRA_China_adj);

% Store results in a table
disp('Business Cycle Statistics (Deseasonalized Data):')
business_cycle_stats_adj = table(["USA"; "France"; "China"], ...
    [volatility_us_adj; volatility_france_adj; volatility_china_adj], ...
    [persistence_us_adj(2); persistence_france_adj(2); persistence_china_adj(2)], ...
    [corr_us_france_adj(1,2); corr_china_france_adj(1,2); corr_china_us_adj(1,2)], ...
    'VariableNames', {'Country', 'Volatility', 'Persistence', 'Correlation'});

disp(business_cycle_stats_adj);

% Plot business cycle of France and China for the whole period
figure;
plot(business_cycle_table_US_France_adj.Quarter, business_cycle_table_US_France_adj.US_Cycle, 'b', 'LineWidth', 1.5); hold on;
plot(business_cycle_table_US_France_adj.Quarter, business_cycle_table_US_France_adj.France_Cycle, 'g', 'LineWidth', 1.5);
xlabel('Year');
ylabel('Business Cycle Component');
title('Quarterly GDP Business Cycle - US & France (Deseasonalized)');
legend('USA', 'France');
grid on;

% Plot before seasonal adjustment
figure;
subplot(2,1,1);
plot(US_France_China_GDP.Quarter, cycle_USA_China, 'b', 'LineWidth', 1.5); hold on;
plot(US_France_China_GDP.Quarter, cycle_FRA_China, 'g', 'LineWidth', 1.5);
plot(US_France_China_GDP.Quarter, cycle_CHN, 'r', 'LineWidth', 1.5);
xlabel('Year');
ylabel('Business Cycle Component');
title('Business Cycle - US, France & China (Before Deseasonalization)');
legend('USA', 'France', 'China');
grid on;

% Plot after seasonal adjustment
subplot(2,1,2);
plot(US_France_China_GDP.Quarter(1:length(cycle_USA_China_adj)), cycle_USA_China_adj, 'b', 'LineWidth', 1.5); hold on;
plot(US_France_China_GDP.Quarter(1:length(cycle_FRA_China_adj)), cycle_FRA_China_adj, 'g', 'LineWidth', 1.5);
plot(US_France_China_GDP.Quarter(1:length(cycle_CHN_adj)), cycle_CHN_adj, 'r', 'LineWidth', 1.5);
xlabel('Year');
ylabel('Business Cycle Component');
title('Business Cycle - US, France & China (After Deseasonalization)');
legend('USA', 'France', 'China');
grid on;

%Volatility Reduction:

%All three countries (USA, France, China) show lower volatility in 
% business cycles after deseasonalization. This suggests that the 
% seasonality in the original data contributed significantly to 
% fluctuations.
%China's volatility remains significantly higher than the USA and 
% France, indicating its economy still experiences greater cyclical 
% swings.

%Increased Persistence:

%The persistence of business cycles (measured by first-order 
% autocorrelation) has increased for all three countries. The USA now 
% has a persistence of 0.912, compared to 0.78 before.
%This suggests that removing seasonal effects has made business cycles 
% appear smoother and more prolonged.

%Stronger Business Cycle Correlations:

%USA-France correlation increased from 0.46985 to 0.49312.
%China-France correlation increased from 0.25548 to 0.43368.
%China-USA correlation increased from 0.17751 to 0.32303.
%These results imply that seasonality was masking some underlying 
% co-movement between economies. Once removed, the true synchronization 
% of business cycles is clearer.

%% 4. Performing Lagged and Lead Correlations for Deseasonalized Data

max_lag = 12; % Define max lag to check (quarters)

% Initialize matrices to store results
corr_lags_adj = zeros(2*max_lag+1, 3); % Rows for -12 to +12 lags

for lag = -max_lag:max_lag
    if lag < 0
        % Leading: Shift second country forward
        corr_lags_adj(lag+max_lag+1, 1) = corr(china_cycle_adj(1:end+lag), cycle_USA_China_adj(-lag+1:end)); 
        corr_lags_adj(lag+max_lag+1, 2) = corr(china_cycle_adj(1:end+lag), cycle_FRA_China_adj(-lag+1:end)); 
        corr_lags_adj(lag+max_lag+1, 3) = corr(us_cycle_adj(1:end+lag), france_cycle_adj(-lag+1:end)); 
    elseif lag > 0
        % Lagging: Shift first country forward
        corr_lags_adj(lag+max_lag+1, 1) = corr(china_cycle_adj(lag+1:end), cycle_USA_China_adj(1:end-lag)); 
        corr_lags_adj(lag+max_lag+1, 2) = corr(china_cycle_adj(lag+1:end), cycle_FRA_China_adj(1:end-lag)); 
        corr_lags_adj(lag+max_lag+1, 3) = corr(us_cycle_adj(lag+1:end), france_cycle_adj(1:end-lag)); 
    else
        % No shift (original correlation)
        corr_lags_adj(lag+max_lag+1, 1) = corr(china_cycle_adj, cycle_USA_China_adj);
        corr_lags_adj(lag+max_lag+1, 2) = corr(china_cycle_adj, cycle_FRA_China_adj);
        corr_lags_adj(lag+max_lag+1, 3) = corr(us_cycle_adj, france_cycle_adj);
    end
end

% Store results in a table
lag_labels = (-max_lag:max_lag)';
lead_lag_table_adj = table(lag_labels, corr_lags_adj(:,1), corr_lags_adj(:,2), corr_lags_adj(:,3), ...
    'VariableNames', {'Lag', 'China_US_Corr', 'China_France_Corr', 'US_France_Corr'});

disp('Lead-Lag Correlations (Deseasonalized Data):')
disp(lead_lag_table_adj)

% Plot Lead-Lag Correlations
figure;
plot(lag_labels, corr_lags_adj(:,1), '-or', 'LineWidth', 1.5); hold on;
plot(lag_labels, corr_lags_adj(:,2), '-ob', 'LineWidth', 1.5);
plot(lag_labels, corr_lags_adj(:,3), '-og', 'LineWidth', 1.5);
xlabel('Lag (quarters)');
ylabel('Correlation');
title('Lead-Lag Correlation of Business Cycles (Deseasonalized Data)');
legend('China-USA', 'China-France', 'USA-France');
grid on;

%Lead-Lag Correlation Improvements:

%The new lead-lag correlation graph looks more structured and smooth 
% compared to the original, which had erratic spikes.

%The curves in the new graph appear more symmetric and display a 
% clearer cyclical pattern, reflecting more consistent relationships 
% between countries’ business cycles.

%The peak of the correlations also appears stronger, reinforcing that 
% deseasonalized data gives a better representation of the actual 
% business cycle comovements.

%% 5. Try the Bry-Boschan (BBQ) Method

% The Bry-Boschan (BBQ) algorithm is a method used to identify turning 
% points (peaks and troughs) in business cycle data. The main goal of the 
% method is to detect economic expansions and recessions by looking at how 
% GDP fluctuates over time.

% The method works by:
% Identifying local maxima as peaks (turning points marking the transition 
% from expansion to recession).
% Identifying local minima as troughs (turning points marking the 
% transition from recession to expansion).
%Ensuring alternation between peaks and troughs (i.e., no two consecutive 
% peaks or troughs).
%Filtering out too-short cycles (e.g., noise) by setting a minimum cycle 
% length (in our case, 8 quarters).

function [peaks, troughs, turning_points] = detect_turning_points(cycle, min_cycle_length, min_cycle_amplitude)
    peaks = [];
    troughs = [];
    
    % Identify initial peaks and troughs
    for i = 2:length(cycle)-1
        if cycle(i) > cycle(i-1) && cycle(i) > cycle(i+1) % Peak condition
            if isempty(peaks) || (i - peaks(end) >= min_cycle_length)
                peaks = [peaks; i];
            end
        elseif cycle(i) < cycle(i-1) && cycle(i) < cycle(i+1) % Trough condition
            if isempty(troughs) || (i - troughs(end) >= min_cycle_length)
                troughs = [troughs; i];
            end
        end
    end

    % Ensure peaks and troughs alternate properly
    valid_peaks = [];
    valid_troughs = [];

    while ~isempty(peaks) && ~isempty(troughs)
        if troughs(1) < peaks(1)
            valid_troughs = [valid_troughs; troughs(1)];
            troughs(1) = []; 
            if ~isempty(peaks)  
                valid_peaks = [valid_peaks; peaks(1)];
                peaks(1) = []; 
            end
        elseif peaks(1) < troughs(1)
            valid_peaks = [valid_peaks; peaks(1)];
            peaks(1) = []; 
            if ~isempty(troughs)  
                valid_troughs = [valid_troughs; troughs(1)];
                troughs(1) = []; 
            end
        else
            peaks(1) = []; 
        end
    end

    % Apply minimum cycle amplitude filtering
    filtered_peaks = [];
    filtered_troughs = [];

    for i = 1:min(length(valid_peaks), length(valid_troughs))
        if abs(cycle(valid_peaks(i)) - cycle(valid_troughs(i))) > min_cycle_amplitude
            filtered_peaks = [filtered_peaks; valid_peaks(i)];
            filtered_troughs = [filtered_troughs; valid_troughs(i)];
        end
    end

    % Assign final alternating peaks & troughs
    peaks = filtered_peaks;
    troughs = filtered_troughs;

    % Ensure the last peak has a corresponding trough
    if ~isempty(peaks) && ~isempty(troughs) && peaks(end) > troughs(end)
        peaks(end) = [];
    end

    % Combine peaks and troughs into a labeled table
    peak_table = table(peaks, repmat("Peak", length(peaks), 1), 'VariableNames', {'Index', 'Type'});
    trough_table = table(troughs, repmat("Trough", length(troughs), 1), 'VariableNames', {'Index', 'Type'});

    % Merge and sort by time
    turning_points = [peak_table; trough_table];
    turning_points = sortrows(turning_points, 'Index');
end

min_cycle_length = 8;   % Minimum 8 quarters between peaks/troughs
min_cycle_amplitude = 0.0075; % Minimum peak-to-trough difference required

[peaks_US, troughs_US, turning_points_US] = detect_turning_points(cycle_USA_adj, min_cycle_length, min_cycle_amplitude);
[peaks_FRA, troughs_FRA, turning_points_FRA] = detect_turning_points(cycle_FRA_adj, min_cycle_length, min_cycle_amplitude);
[peaks_CHN, troughs_CHN, turning_points_CHN] = detect_turning_points(cycle_CHN_adj, min_cycle_length, min_cycle_amplitude);

% Display the Corrected Peaks and Troughs
disp('Peaks & Troughs:');
disp('US Peaks:'); disp(peaks_US);
disp('US Troughs:'); disp(troughs_US);
disp('France Peaks:'); disp(peaks_FRA);
disp('France Troughs:'); disp(troughs_FRA);
disp('China Peaks:'); disp(peaks_CHN);
disp('China Troughs:'); disp(troughs_CHN);

%% 6. Compute Business Cycle Statistics from BBQ Method

function [expansion_duration, recession_duration, mean_exp_growth, mean_rec_growth] = compute_durations(peaks, troughs, cycle)
    % Ensure at least one peak and one trough exist
    if isempty(peaks) || isempty(troughs)
        expansion_duration = NaN;
        recession_duration = NaN;
        mean_exp_growth = NaN;
        mean_rec_growth = NaN;
        return;
    end

% The function compute_durations takes three inputs:
% peaks: The list of peak indices (economic highs).
% troughs: The list of trough indices (economic lows).
% cycle: The business cycle fluctuations (filtered GDP series).
% The function returns:
% expansion_duration: The number of periods (quarters) between each peak 
% and the next trough (length of economic expansions).
% recession_duration: The number of periods between each trough and the 
% next peak (length of economic recessions).
% mean_exp_growth: The average GDP growth during expansion phases.
% mean_rec_growth: The average GDP decline during recessions.

% The first check ensures that we have at least one peak and one trough 
% in the detected turning points.
%If not, the function returns NaN values (meaning we don’t have enough 
% data to compute cycle durations).

    % Ensure alternation: First peak should be followed by a trough
    if peaks(1) < troughs(1)
        expansion_duration = diff(peaks); % Expansion duration = time between peaks
        recession_duration = diff(troughs); % Recession duration = time between troughs
    else
        recession_duration = diff(troughs);
        expansion_duration = diff(peaks);
    end

% If the first turning point is a peak, we compute:
% Expansion durations from peak to peak.
% Recession durations from trough to trough.
% If the first turning point is a trough, the assignments are reversed.

    % Compute mean expansion and recession growth rates
    mean_exp_growth = mean(cycle(peaks));
    mean_rec_growth = mean(cycle(troughs));
end

% mean_exp_growth: The average GDP cycle value at peaks (expansion).
% mean_rec_growth: The average GDP cycle value at troughs (recession).

% Compute statistics for each country
[duration_expansion_US, duration_recession_US, mean_expansion_US, mean_recession_US] = ...
    compute_durations(peaks_US, troughs_US, cycle_USA_adj);

[duration_expansion_FRA, duration_recession_FRA, mean_expansion_FRA, mean_recession_FRA] = ...
    compute_durations(peaks_FRA, troughs_FRA, cycle_FRA_adj);

[duration_expansion_CHN, duration_recession_CHN, mean_expansion_CHN, mean_recession_CHN] = ...
    compute_durations(peaks_CHN, troughs_CHN, cycle_CHN_adj);

% Now that we have the compute_durations function, we apply it to 
% the USA, France, and China.
% Each call to compute_durations gives us:
% How long expansions and recessions last.
% How strong expansions and recessions are.

% Compute average durations
duration_expansion_US_mean = mean(duration_expansion_US);
duration_recession_US_mean = mean(duration_recession_US);

duration_expansion_FRA_mean = mean(duration_expansion_FRA);
duration_recession_FRA_mean = mean(duration_recession_FRA);

duration_expansion_CHN_mean = mean(duration_expansion_CHN);
duration_recession_CHN_mean = mean(duration_recession_CHN);

% We compute the mean duration of expansions and recessions for each
% country

% Store business cycle statistics in a table
business_cycle_durations = table(["USA"; "France"; "China"], ...
    [duration_expansion_US_mean; duration_expansion_FRA_mean; duration_expansion_CHN_mean], ...
    [duration_recession_US_mean; duration_recession_FRA_mean; duration_recession_CHN_mean], ...
    [mean_expansion_US; mean_expansion_FRA; mean_expansion_CHN], ...
    [mean_recession_US; mean_recession_FRA; mean_recession_CHN], ...
    'VariableNames', {'Country', 'Avg_Expansion_Duration', 'Avg_Recession_Duration', ...
                      'Mean_Expansion_Growth', 'Mean_Recession_Growth'});

% Display results
disp('Business Cycle Durations:');
disp(business_cycle_durations);

% This table summarizes:
%The average number of quarters in expansions and recessions.
%How much GDP grows or contracts during these phases.

%% 7. Plot Business Cycles for France & USA (Full Period)

figure;
hold on;

% Adjust x-axis to match available cycle lengths
plot(US_France_GDP.Quarter(end-length(cycle_USA_adj)+1:end), cycle_USA_adj, 'b', 'LineWidth', 1.5); % USA
plot(US_France_GDP.Quarter(end-length(cycle_FRA_adj)+1:end), cycle_FRA_adj, 'g', 'LineWidth', 1.5); % France

% Mark Turning Points (Peaks & Troughs)
scatter(US_France_GDP.Quarter(peaks_US), cycle_USA_adj(peaks_US), 50, 'bo', 'filled'); % USA Peaks
scatter(US_France_GDP.Quarter(troughs_US), cycle_USA_adj(troughs_US), 50, 'bs', 'filled'); % USA Troughs

scatter(US_France_GDP.Quarter(peaks_FRA), cycle_FRA_adj(peaks_FRA), 50, 'go', 'filled'); % France Peaks
scatter(US_France_GDP.Quarter(troughs_FRA), cycle_FRA_adj(troughs_FRA), 50, 'gs', 'filled'); % France Troughs

% Formatting
xlabel('Year');
ylabel('Business Cycle Component');
title('Business Cycles - USA & France (Full Period)');
legend('USA Cycle', 'France Cycle', ...
       'USA Peaks', 'USA Troughs', 'France Peaks', 'France Troughs');
grid on;
hold off;

%% 8. Plot Business Cycles for France, USA, and China (Full Period)

figure;
hold on;

% Adjust x-axis to match available cycle lengths
plot(US_France_GDP.Quarter(end-length(cycle_USA_adj)+1:end), cycle_USA_adj, 'b', 'LineWidth', 1.5); % USA
plot(US_France_GDP.Quarter(end-length(cycle_FRA_adj)+1:end), cycle_FRA_adj, 'g', 'LineWidth', 1.5); % France
plot(US_France_China_GDP.Quarter(end-length(cycle_CHN_adj)+1:end), cycle_CHN_adj, 'r', 'LineWidth', 1.5); % China

% Mark Turning Points (Peaks & Troughs)
scatter(US_France_GDP.Quarter(peaks_US), cycle_USA_adj(peaks_US), 50, 'bo', 'filled'); % USA Peaks
scatter(US_France_GDP.Quarter(troughs_US), cycle_USA_adj(troughs_US), 50, 'bs', 'filled'); % USA Troughs

scatter(US_France_GDP.Quarter(peaks_FRA), cycle_FRA_adj(peaks_FRA), 50, 'go', 'filled'); % France Peaks
scatter(US_France_GDP.Quarter(troughs_FRA), cycle_FRA_adj(troughs_FRA), 50, 'gs', 'filled'); % France Troughs

scatter(US_France_China_GDP.Quarter(peaks_CHN), cycle_CHN_adj(peaks_CHN), 50, 'ro', 'filled'); % China Peaks
scatter(US_France_China_GDP.Quarter(troughs_CHN), cycle_CHN_adj(troughs_CHN), 50, 'rs', 'filled'); % China Troughs

% Formatting
xlabel('Year');
ylabel('Business Cycle Component');
title('Business Cycles - USA, France & China (Full Period)');
legend('USA Cycle', 'France Cycle', 'China Cycle', ...
       'USA Peaks', 'USA Troughs', 'France Peaks', 'France Troughs', 'China Peaks', 'China Troughs');
grid on;
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% ANNUAL DATA ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 1. Loading the dataset

% Setting the variable names of the data set
opts = detectImportOptions('Penn_World_Tables.csv'); % Auto-detect import settings
opts.DataLines = 2;  % Ensure data starts from row 2
opts.VariableNames = ["ISO_code", "Country", "Variable_code", "Variable_name", ...
    "1950", "1951", "1952", "1953", "1954", "1955", "1956", "1957", "1958", "1959", ...
    "1960", "1961", "1962", "1963", "1964", "1965", "1966", "1967", "1968", "1969", ...
    "1970", "1971", "1972", "1973", "1974", "1975", "1976", "1977", "1978", "1979", ...
    "1980", "1981", "1982", "1983", "1984", "1985", "1986", "1987", "1988", "1989", ...
    "1990", "1991", "1992", "1993", "1994", "1995", "1996", "1997", "1998", "1999", ...
    "2000", "2001", "2002", "2003", "2004", "2005", "2006", "2007", "2008", "2009", ...
    "2010", "2011", "2012", "2013", "2014", "2015", "2016", "2017", "2018", "2019"]; % Manually setting variable names

data_table = readtable('Penn_World_Tables.csv', opts); % Read table with corrected settings

% Display the first few rows to verify
disp(data_table(1:5,1:6))

%% 2. Data cleaning

% Count missing values per column
missing_counts = sum(ismissing(data_table));

% Display missing values count for each column
disp(table(data_table.Properties.VariableNames', missing_counts', 'VariableNames', {'Column', 'MissingValues'}))

% Select only years from 1970 onward
year_columns = 5 + (1970-1950):width(data_table); % Identify columns from 1970 onwards
data_table = data_table(:, [1:4, year_columns]); % Keep categorical variables and selected years

% Display first few rows to check
disp(data_table(1:5,1:6))

%% 3. Validating data types

% Check data types of each column
disp(varfun(@class, data_table))

% Converting cell arrays to categorical variables
data_table.ISO_code = categorical(data_table.ISO_code);
data_table.Country = categorical(data_table.Country);
data_table.Variable_code = categorical(data_table.Variable_code);
data_table.Variable_name = categorical(data_table.Variable_name);

% Check data types of each column
disp(varfun(@class, data_table))

%% 4. Extracting and Structuring GDP Data

% Define Countries and GDP Variable
selected_countries = categorical(["CHN", "USA", "FRA"]); % Select China, US, and France
gdp_variable = categorical("cgdpo"); % Output-side Real GDP (PPP-adjusted)

%Extract GDP Data for Selected Countries
gdp_data = data_table(ismember(data_table.ISO_code, selected_countries) & ...
                      data_table.Variable_code == gdp_variable, :);

% Display the extracted data to verify
disp(gdp_data)

% Convert GDP Table to Matrix
% Extract only numeric GDP data (years 1970-2019)
gdp_matrix = table2array(gdp_data(:, 5:end)); % Columns 5:end are years

% Ensure GDP matrix is transposed (years as rows, countries as columns)
gdp_matrix = gdp_matrix';

% Define correct year labels as a column vector
years = (1970:2019)'; 

% Ensure correct display of years and GDP matrix
gdp_table = table(years, gdp_matrix(:,1), gdp_matrix(:,2), gdp_matrix(:,3), ...
    'VariableNames', {'Year', 'China_GDP', 'US_GDP', 'France_GDP'});

% Display the table properly
disp(gdp_table)

%% 5. Convert GDP to Logarithms

% Convert GDP values to log format
log_gdp = log(gdp_matrix);

% Store log-transformed GDP in a structured table for reference
log_gdp_table = table(years, log_gdp(:,1), log_gdp(:,2), log_gdp(:,3), ...
    'VariableNames', {'Year', 'Log_China_GDP', 'Log_US_GDP', 'Log_France_GDP'});

% Display first 5 rows to check
disp(log_gdp_table)

%% 6. Apply the HP Filter

lambda = 100; % Standard smoothing parameter for annual data
% (as per Ravn & Uhlig, 2002)

% Apply the HP filter to each country
[cycle_CHN, trend_CHN] = hpfilter(log_gdp(:,1), lambda);
[cycle_USA, trend_USA] = hpfilter(log_gdp(:,2), lambda);
[cycle_FRA, trend_FRA] = hpfilter(log_gdp(:,3), lambda);

% Store results in a single table
business_cycle_table = table(years, cycle_CHN, cycle_USA, cycle_FRA, ...
    'VariableNames', {'Year', 'China_Cycle', 'US_Cycle', 'France_Cycle'});

% Display first few rows to verify
disp(business_cycle_table)

%% 7. Plot the Business Cycle

% Plot the business cycle
figure;
plot(years, business_cycle_table.China_Cycle, 'r', 'LineWidth', 1.5); hold on;
plot(years, business_cycle_table.US_Cycle, 'b', 'LineWidth', 1.5);
plot(years, business_cycle_table.France_Cycle, 'g', 'LineWidth', 1.5);
xlabel('Year');
ylabel('Business Cycle Component');
title('GDP Business Cycle Fluctuations');
legend('China', 'USA', 'France');
grid on;

%% 8. Compute Business Cycle Statistics

% Compute Business Cycle Statistics

% Extract business cycle components
china_cycle = business_cycle_table.China_Cycle;
us_cycle = business_cycle_table.US_Cycle;
france_cycle = business_cycle_table.France_Cycle;

% Compute Volatility (Standard Deviation)
volatility_china = std(china_cycle);
volatility_us = std(us_cycle);
volatility_france = std(france_cycle);

% Compute Persistence (First-Order Autocorrelation)
persistence_china = autocorr(china_cycle, 1);
persistence_us = autocorr(us_cycle, 1);
persistence_france = autocorr(france_cycle, 1);

% Compute Co-movement (Correlation between Countries)
corr_china_us = corrcoef(china_cycle, us_cycle);
corr_china_france = corrcoef(china_cycle, france_cycle);
corr_us_france = corrcoef(us_cycle, france_cycle);

% Display Results
disp('Business Cycle Statistics:')
business_cycle_stats = table(["China"; "USA"; "France"], ...
    [volatility_china; volatility_us; volatility_france], ...
    [persistence_china(2); persistence_us(2); persistence_france(2)], ...
    [corr_china_us(1,2); corr_china_france(1,2); corr_us_france(1,2)], ...
    'VariableNames', {'Country', 'Volatility', 'Persistence', 'Correlation'});

disp(business_cycle_stats)

%Volatility:
% China has the highest volatility (0.0649). This means China's GDP
% fluctuations are the most pronounced.
% USA is in the middle (0.0436). The US has moderate business cycle
% fluctuations.
%France has the lowest volatility (0.0250) → France's GDP fluctuations
% are the smallest, indicating more stability.
% China's economy experiences larger boom and bust cycles, while France's
% economy appears to be more stable.

% Persistence (First-Order Autocorrelation):
%USA has the highest persistence (0.7969). Meaning that its business cycle
% deviations tend to be more prolonged over time.
%China (0.7579) and France (0.7341) have lower persistence but still show
% some degree of continuity.
%If there is a recession or an expansion, it is likely to last longer in
% the USA than in China or France.

%Co-movement (Correlation):
%France has the strongest correlation (0.4534) with the other countries.
% Suggesting France's economy is more synchronized with global economic trends.
%China (-0.1907) and the USA (-0.0395) show weak or negative correlations.
% Meaning their business cycles do not move closely together.
%The weak or even negative correlation between China and the USA may
% indicate that their business cycles are influenced by different economic
% forces (e.g., policy, trade, global shocks).
%France being more positively correlated suggests it is more integrated
% with the global economic system.

%% 9. Double Checking Correlations with Lagged and Lead Correlations

% Check Lead-Lag Correlations Between Countries

max_lag = 3; % Define max lag to check

% Initialize matrices to store results
corr_lags = zeros(2*max_lag+1, 3); % Rows for -3 to +3 lags, columns for
% each country pair

for lag = -max_lag:max_lag
    if lag < 0
        % Leading: Shift second country forward
        corr_lags(lag+max_lag+1, 1) = corr(china_cycle(1:end+lag), us_cycle(-lag+1:end)); 
        corr_lags(lag+max_lag+1, 2) = corr(china_cycle(1:end+lag), france_cycle(-lag+1:end)); 
        corr_lags(lag+max_lag+1, 3) = corr(us_cycle(1:end+lag), france_cycle(-lag+1:end)); 
    elseif lag > 0
        % Lagging: Shift first country forward
        corr_lags(lag+max_lag+1, 1) = corr(china_cycle(lag+1:end), us_cycle(1:end-lag)); 
        corr_lags(lag+max_lag+1, 2) = corr(china_cycle(lag+1:end), france_cycle(1:end-lag)); 
        corr_lags(lag+max_lag+1, 3) = corr(us_cycle(lag+1:end), france_cycle(1:end-lag)); 
    else
        % No shift (original correlation)
        corr_lags(lag+max_lag+1, 1) = corr(china_cycle, us_cycle);
        corr_lags(lag+max_lag+1, 2) = corr(china_cycle, france_cycle);
        corr_lags(lag+max_lag+1, 3) = corr(us_cycle, france_cycle);
    end
end

% Store results in a table
lag_labels = (-max_lag:max_lag)';
lead_lag_table = table(lag_labels, corr_lags(:,1), corr_lags(:,2), corr_lags(:,3), ...
    'VariableNames', {'Lag', 'China_US_Corr', 'China_France_Corr', 'US_France_Corr'});

% Display lead-lag correlations
disp('Lead-Lag Correlations:')
disp(lead_lag_table)

% Plot Lead-Lag Correlations
figure;
plot(lag_labels, corr_lags(:,1), '-or', 'LineWidth', 1.5); hold on;
plot(lag_labels, corr_lags(:,2), '-ob', 'LineWidth', 1.5);
plot(lag_labels, corr_lags(:,3), '-og', 'LineWidth', 1.5);
xlabel('Lag (years)');
ylabel('Correlation');
title('Lead-Lag Correlation of Business Cycles');
legend('China-USA', 'China-France', 'USA-France');
grid on;

%% 10. Notes on Historical Events Affecting Business Cycles

%_____________________________________________________________________%
% China: High Volatility (0.0649), Moderate Persistence(0.7579), Weak
% Correlation (-0.1907) with the USA
%_____________________________________________________________________%

% 1978–1990: Economic Reforms and Opening Up
% After Deng Xiaoping’s reforms in 1978, China shifted from a centrally
% planned economy to market-oriented reforms, causing high GDP growth.The
% shift resulted in large fluctuations in economic activity due to rapid
% industrialization.

%1990s–2000s: Export-Led Growth and WTO Accession (2001)
%China’s economic boom was fueled by exports and foreign direct investment
% (FDI). Business cycle volatility increased as the economy transitioned
% from state-led to private-sector-driven growth.

%2008–2009: Global Financial Crisis
%China launched a massive stimulus package (~$586 billion) to counter the
% global crisis. Business cycle fluctuations became less synchronized with
% Western economies due to state intervention.

%2015–2016: Stock Market Crash and Slowdown
%China faced a stock market crash and concerns over a slowing economy.This 
% period reinforced volatility, as the government shifted between market 
% liberalization and intervention.

%China’s high business cycle volatility comes from policy-driven growth, 
% rapid transitions, and external trade shocks. Its weak correlation with 
% the USA suggests China’s cycles are more influenced by domestic policies
% than global trends.

%________________________________________________________________________%
% USA: Moderate Volatility(0.0436), High Persistence (0.7969 so recessions
% and expansions last longer), Weak Correlation (-0.0395) with China and
% France.
%________________________________________________________________________%

%1970s: Stagflation and Oil Crises (1973, 1979)
%The US economy suffered from high inflation and stagnation, leading to 
% prolonged slowdowns.Persistence was high, as economic downturns lasted 
%for several years.

%1980s: Reaganomics and Recovery
%Tax cuts, deregulation, and monetary policy led to a long expansion.
% Business cycles had high persistence—once growth started, it lasted for
% a decade.

%2000s: Dot-Com Bubble and 2008 Financial Crisis
%The 2000 dot-com crash and 2008 housing bubble burst led to recessions 
% with prolonged recoveries. The Great Recession (2008-2009) saw the 
% largest contraction in US GDP since the Great Depression.

%2020: COVID-19 Pandemic
%The sharpest GDP drop in history followed by a strong rebound. 
% Unprecedented fiscal and monetary policy responses led to high
% volatility.

%The high persistence of the US business cycle is due to financial
% market-driven expansions and deep recessions. The weak correlation
% with China and France suggests the US follows its own financial-driven
% cycle rather than external trade shocks.

%__________________________________________________________________%
% France: Lowest Volatility(0.0250), Moderate Persistence (0.7341),
% Stronger Correlation (0.4534) with global business cycles.
%__________________________________________________________________%

%1970s: Oil Crises and European Recession
%France suffered from rising oil prices and slow industrial growth. This
% period marked the beginning of lower volatility, as the government
% focused on stability over rapid expansion.

%1992–1993: European Exchange Rate Mechanism (ERM) Crisis
%The European monetary system collapse led to economic instability. France
% became more aligned with the European business cycle.

%2000s: Euro Adoption and Financial Crises
%France’s adoption of the euro (1999) integrated it further into the
% Eurozone economy. 2008 Financial Crisis and Eurozone Crisis (2010-2012)
% caused GDP contractions, but France remained less volatile than the USA
% due to strong social policies.

%2019–2022: COVID-19 and Inflation Crisis
%The French economy contracted less sharply than the USA due to government
% intervention (subsidies, furlough schemes). Business cycle correlation
% with global trends increased due to Eurozone dependency.

%France’s low volatility results from strong state intervention and
% Eurozone policies. The higher correlation (0.4534) with the US and China
% suggests that France follows global economic trends more than China.

