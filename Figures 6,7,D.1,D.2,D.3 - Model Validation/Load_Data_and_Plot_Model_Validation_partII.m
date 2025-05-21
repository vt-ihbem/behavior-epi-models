%% Load and Plot Data for Model Validation (part II)
clear; close all; clc

%% Load Data and Simulated Data
Data = readtable('Compaing 6 models.xlsx');

%% Plot Data - Appendix Figure D2

names = Data.Properties.VariableNames;
names{6} = 'SEIRS + seasonality';
names{8} = 'SEIRb + seasonality';
names{9} = 'SEIRSb + seasonality';

style = {'-',':','--','-.','-.','--','-',};
f1 = figure(1);
order = [3 9 8 7 6 5 4];
for j = order
    plot(Data.Date,Data.(j),'linewidth',2,'linestyle',style{j-2})
    hold on
end
hold off
l = legend(names{order});
set(l,'box','off','location','northwest')
set(gca,'fontsize',18)

ylabel('Daily death (7 day average)','fontsize',18)

%% Save plot

return
saveas(f1,'FigureD2.png')