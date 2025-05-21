%% Load and Plot Data for Model Validation (part I)
clear; close all; clc

%% Load Data and Simulated Data

Data = readtable('Model Validation Data.xlsx','Sheet','Data (nominal)');
Sim = readtable('Model Validation Data.xlsx','Sheet','Simulation (nominal)');

%% Plot Data - Appendix Figure D.3

names = Data.Properties.VariableNames;
names{11} = 'District Of Columbia';
names{32} = 'New Hampshire';
names{33} = 'New Jersey';
names{34} = 'New Mexico';
names{35} = 'New York';
names{36} = 'North Carolina';
names{37} = 'North Dakota';
names{42} = 'Rhode Island';
names{43} = 'South Carolina';
names{44} = 'South Dakota';
names{51} = 'West Virginia';

fD3 = figure(2);
set(fD3,'Position',[1 1 1200 1000])
clf
t = tiledlayout(9,6);
for j = 3:54
    nexttile
    plot(Data.Date,Data.(j),'linewidth',1.5)
    hold on
    plot(Sim.Date,Sim.(j),'linewidth',2)
    hold off
    set(gca,'fontsize',12)
    if j<49
        set(gca,'xticklabel',[])
    end
    title(names{j})
end
ylabel(t,'Daily death (7 day average)','fontsize',16)

%% Plot Data - Figure 7

id(1) = find(strcmp(names,'New York'));
id(2) = find(strcmp(names,'Texas'));
id(3) = find(strcmp(names,'California'));
id(4) = find(strcmp(names,'Hawaii'));
id(5) = find(strcmp(names,'Utah'));
id(6) = find(strcmp(names,'South Dakota'));
id(7) = find(strcmp(names,'Vermont'));
id(10) = find(strcmp(names,'Wyoming'));
id(8) = find(strcmp(names,'USA'));

f7 = figure(7);
clf
set(f7,'Position',[1 1 800 1000])
for j = [1:8 10]
    subplot(4,3,j)

    if j==8
        subplot(4,3,[8:9 11:12])
    end
    plot(Data.Date,Data.(id(j)),'linewidth',1.5)
    hold on
    plot(Sim.Date,Sim.(id(j)),'linewidth',2)
    hold off
    if j <=7
        set(gca,'xticklabel','')
    end
    if j == 8
        l = legend('Data','Model');
        set(l,'box','off','location','northwest')
    end

    set(gca,'fontsize',16)
    title(names{id(j)})
end
han=axes(f7,'visible','off');
han.YLabel.Visible='on';
ylabel(han,'Daily death (7 day average)','fontsize',16);


%% Load Data - R-squared and Absolute Mean error

SummaryData = readtable('Model Validation Summary.xlsx');
Rsq = SummaryData{1:51,3:7};
MSE = SummaryData{1:51,8:12};
Names = SummaryData{1:51,2};

%% Plot Data - Figure 6 

types = {'SEIRS','SEIRS with seasonality','SEIRSb','SEIRb with seasonality','SEIRSb with seasonality'};
f6 = figure(6);
boxplot(Rsq,'Whisker',1.5,'Widths',0.5,'symbol','')
ylim([0 1])
set(gca,'xticklabel',types)
ylabel('R squared')
set(gca,'fontsize',18)

%% Plot Data - Appendix Figure 

fD1 = figure(4);
boxplot(MSE,'Whisker',1.5,'symbol','') 
set(gca,'xticklabel',types)
ylabel('Mean absolute error')
set(gca,'fontsize',18)
ylim([0 6.5])
h=findobj('LineStyle','--'); set(h, 'LineStyle','-');

%% Save plots

return
saveas(f6,'Figure6.png')
saveas(fD1,'FigureD1.png')
saveas(fD3,'FigureD3.png')
saveas(f7,'Figure7.png')

%% To get sheet names automatically
filename = 'Model Validation Data.xlsx';
sheetNames = sheetnames(filename);