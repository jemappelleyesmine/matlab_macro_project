%%%%%%%%%%%%%%%%%%%%%%%%%%
% MACROECONOMICS PROJECT %
%%%%%%%%%%%%%%%%%%%%%%%%%%


cd '/Users/yesminehachana/Desktop/Classes/Dauphine/2nd Semester/Macroeconomics MATLAB Project/Repository'

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

% Comments on business cycle stats

%Volatility
% China has the highest volatility (0.0649). This means China's GDP
% fluctuations are the most pronounced.
% USA is in the middle (0.0436). The US has moderate business cycle
% fluctuations.
%France has the lowest volatility (0.0250) → France's GDP fluctuations
% are the smallest, indicating more stability.
% China's economy experiences larger boom and bust cycles, while France's
% economy appears to be more stable.

% Persistence (First-Order Autocorrelation)
%USA has the highest persistence (0.7969). Meaning that its business cycle
% deviations tend to be more prolonged over time.
%China (0.7579) and France (0.7341) have lower persistence but still show
% some degree of continuity.
%If there is a recession or an expansion, it is likely to last longer in
% the USA than in China or France.

%Co-movement (Correlation)
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

%_____________________________________________________________________
% China: High Volatility (0.0649), Moderate Persistence(0.7579), Weak
% Correlation (-0.1907) with the USA
%_____________________________________________________________________

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

%_________________________________________________________________________
% USA: Moderate Volatility(0.0436), High Persistence (0.7969 so recessions
% and expansions last longer), Weak Correlation (-0.0395) with China and
% France.
%_________________________________________________________________________

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

%__________________________________________________________________
% France: Lowest Volatility(0.0250), Moderate Persistence (0.7341),
% Stronger Correlation (0.4534) with global business cycles.
%__________________________________________________________________

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


