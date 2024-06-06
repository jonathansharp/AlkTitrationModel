%% Example TITRATE.m function run
%
% This example function uses TITRATE.m (which is downloaded from
% (github.com/jonathansharp/AlkTitrationModel) to simulate total alkalinity
% (TA) titrations on a profile of open-ocean samples. We will assume a
% uniform total organic acid concentration of 9 umol/kg, composed of three
% separate acids, each with concentrations of 3 umol/kg and pKs equal to 4,
% 5, and 6. We will then plot the difference between modelled TA values
% determined using different titration methods and inorganic alkalinity.
% These differences will represent the excess alkalinity amount detected by
% each method.

% This is data from Station 156 of the P16N cruise in 2015 (44' N, 155' W)
PRS    = [3.4 40.6 69.5 109.1 174.6 250.0 332.7 432.4 550.1 701.0 850.9 1000.2 1167.2 1366.8 1599.9 1933.4 2366.6 2866.0 3366.6 3865.3 4367.8 4866.4 5282.4 5727.4];
TMP    = [11.805 10.0835 9.0949 8.7028 8.1385 7.2186 6.0406 4.8193 4.2995 3.9305 3.4850 3.1329 2.8238 2.5597 2.3071 2.0166 1.7961 1.6101 1.5137 1.4757 1.4998 1.5522 1.6036 1.6633];
SAL    = [33.095 33.0915 33.1315 33.7722 33.9339 33.942 33.9295 33.9477 34.053 34.1861 34.2706 34.3405 34.4022 34.4614 34.5131 34.5741 34.6161 34.6497 34.6675 34.6809 34.6850 34.6858 34.6859 34.6873];
SIL    = [4.21 7.44 8.51 15.06 25.43 38.92 53.79 72.66 89.96 106.57 121.23 132.96 142.73 151.81 159.73 165.97 167.73 164.21 159.80 157.56 160.39 162.83 163.42 163.52];
PO4    = [0.61 0.74 0.84 1.05 1.38 1.76 2.16 2.57 2.82 2.97 3.06 3.12 3.11 3.11 3.08 3.00 2.87 2.73 2.62 2.55 2.52 2.52 2.52 2.51];
DIC    = [2021.9 2028.9 2047.0 2096.7 2135.5 2179.2 2223.2 2269.0 2306.4 2338.4 2360.0 2375.0 2384.2 2393.9 2401.2 2387.7 2382.0 2362.0 2348.2 2338.6 2337.5 2337.1 2337.2 2338.8];
PH     = [7.8744 7.8508 7.8126 7.744 7.6856 7.5933 7.4959 7.3948 7.3341 7.3047 7.2910 7.2882 7.2968 7.2983 7.3124 7.3600 7.4165 7.4772 7.5225 7.5522 7.5653 7.5687 7.5691 7.5689];
PH_TMP = [25.00 24.98 25.00 24.99 25.00 25.01 25.00 25.00 24.98 25.00 25.00 25.00 24.97 25.01 25.00 24.99 24.99 24.99 24.99 25.00 24.98 24.98 24.99 25.01];
TA_M   = [2216.6 2215.2 2216.7 2250.9 2261.1 2278.2 2280.7 2301.5 2320.7 2336.8 2355.3 2369.7 2381.0 2393.9 2405.2 2413.0 2420.8 2422.2 2423.4 2425.4 2429.8 2431.6 2431.0 2432.4];
% Make CO2 system calculations
[PH_20,headers]  = CO2SYS(PH,DIC,3,2,SAL,PH_TMP,20,0,0,SIL,PO4,0,0,1,10,1,2,2);
idx = find(strcmp(headers,'pHout')); PH_20 = PH_20(:,idx);

% This defines three separate organic acids with concentrations of 3
% umol/kg each and pKs of 4, 5, and 6
TORG1  = 3;
TORG2  = 3;
TORG3  = 3;
pKOrg1 = 4;
pKOrg2 = 5;
pKOrg3 = 6;

% This runs the titration model for each sample, assuming a sample mass of
% 200 grams, an acid concentration of 0.1 mol/kg, and a titration
% temperature of 20 C
[Result,Result_Headers] = TITRATE(SAL,20,DIC,PH_20,PO4,SIL,0,0,200,0.1,...
                             TORG1,pKOrg1,TORG2,pKOrg2,TORG3,pKOrg3);

% This plots excess alkalinity measured by each method as a function of
% depth. Markers are color-coded by measured alkalinity.
figure; clf; hold on
scatter(Result(:,6)-Result(:,1),PRS,70,TA_M,'Filled','Marker','^');
scatter(Result(:,7)-Result(:,1),PRS,70,TA_M,'Filled','Marker','s');
scatter(Result(:,4)-Result(:,1),PRS,70,TA_M,'Filled','Marker','o');
scatter(Result(:,3)-Result(:,1),PRS,70,TA_M,'Filled','Marker','d');
xlim([3 9]); % set x-axis limits
set(gcf,'Position',[200,200,600,600]);
set(gca,'Ydir','reverse'); % reverse y-axis direction
set(gca,'Fontsize',15);
xlabel('TA_{model} - TA_{inorg} (\mumol/kg)');
ylabel('Depth (dbar)');
legend({'Diff. Deriv.','Single-Step','Closed-Cell','Open-Cell'},...
    'Location','northoutside','NumColumns',4); % add legend
c = colorbar('northoutside'); c.Label.String = 'TA_{meas} (\mumol/kg)';
hold off