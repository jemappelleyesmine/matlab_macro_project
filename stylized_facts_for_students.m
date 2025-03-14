
% Getting Business Cycle stylized facts 
% OECD Stats as Source
% quarterly real output

% --------------------------------------------
% --- Load the data

clear all
close all
clc

%%%% GDP series

% Read the data: GDP series
% Already countries in column / quarter by line

%database = transpose(database);
%country_list= transpose(txt(8:17,1));
%period = transpose(txt(6,4:end));

[database, txt]  = xlsread('OECD_MEI_Synthesis_v2.xlsx','GDP'); 
country_list= txt(7, 2:16);
period = txt(8:end, 1);
% From 1970:1 to 2019:3 for all countries 





% Take the series in log
database(:,:)=log(database(:,:));
lambda=1600;

% Nb of countries, periods, matrix definition
nbc  = size(database,2); % nb of countries
nbp  = size(database,1); % nb of periods
cycle = zeros(nbp,nbc);
trend = zeros(nbp,nbc);

% time period
tq=[1970:1/4:2018];

for jj=1:nbc;
    [cycle(:,jj), trend(:,jj)]=hpfilter(database(:,jj),lambda); 
    
end;



% Ordre des pays
% 1 'AUS' 2 'AUT' 3 'BEL' 4 'CAN' 5 'FRA' 6 'DEU' 7 'IRL' 8 'ITA' 9 'JPN' 10 'KOR' 11 'NLD' 12 'PRT' 13 'ESP' 14 'GBR' 15 'ZAF'


% Cycle, for a sample of countries
figure;
plot(tq,cycle(:,8),'r'); hold on;plot(tq,cycle(:,7),'g');hold on; plot(tq,cycle(:,5),'b');hold on; plot(tq,cycle(:,14),'k');hold on; plot(tq,cycle(:,12),'y');hold on; plot(tq,cycle(:,6),'b');
title('GDP Business cycle component');hold on; plot(tq,cycle(:,6),'b'); 
xlabel('quarters');
legend('ITA','IRL','FRA','GBR', 'PRT', 'DEU');



% Cycle - Trend for one country, say FRANCE
ref=zeros(nbp,1);
figure;
%subplot(2,1,1);
plot(tq,cycle(:,7),'r');hold on;plot(tq,ref,':')
xlabel('quarter'); legend('France, Business Cycle');

figure;
%subplot(2,1,2);
plot(tq,trend(:,7),'r');
xlabel('quarter'); legend('France, Trend');



%%%% Stylized facts for a list of countries

filename = 'stylized_facts_largesample.xlsx';
% the excel file to export the results


% Set of countries selected, say France and Germany
country_select=strvcat('FRA','DEU');

data_select=zeros(nbp,6); % nb of macro variables in column
% say, you are interested in 6 series

for ii = 1:size(country_select,1);

    if ii==1; country_name='FRA'; jj=5; end;
    if ii==2; country_name='DEU'; jj=6; end;
   
 
    %%% TO BE DONE %%%
    
    % load the appropriate data series
    
    % extract the HP component
      
    % plot series if you want
    
    % determine the second-order moments

    % Hint: Some useful functions for that
    % sig = 100*std(cycle); For volatility

    % Instantaneous correlation
    %correl0 = corrcoef(cycle);

    % To calculate autocorrelation, correlation with lead-lag
    % one way is to build the serie in t that corresponds to the values in
    % t-1, then apply corrcoef
    
    
    
end;





