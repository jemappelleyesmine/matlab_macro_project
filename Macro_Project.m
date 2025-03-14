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

lambda = 1600; % Standard smoothing parameter for quarterly data

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

