clear all; close all; clc
%% Plots: PRCC and p-values

load('LHS_PRCC_SEIRS_output.mat','sortedPRCC_1','sortedPRCC_2','sortedPRCC_3','sortedstat_p1','sortedstat_p2','sortedstat_p3','sortedparams','pars')
SsortedPRCC_1 = sortedPRCC_1;
SsortedPRCC_2 = sortedPRCC_2;
SsortedPRCC_3 = sortedPRCC_3;
Ssortedstat_p1 = sortedstat_p1; 
Ssortedstat_p2 = sortedstat_p2; 
Ssortedstat_p3 = sortedstat_p3; 
Ssortedparams = sortedparams;

load('LHS_PRCC_SEIR_output.mat','sortedPRCC_1','sortedPRCC_2','sortedPRCC_3','sortedstat_p1','sortedstat_p2','sortedstat_p3')
sortedPRCC_1 = [sortedPRCC_1 0];
sortedPRCC_2 = [sortedPRCC_2 0];
sortedPRCC_3 = [sortedPRCC_3 0];
sortedstat_p1 = [sortedstat_p1 1]; 
sortedstat_p2 = [sortedstat_p2 1]; 
sortedstat_p3 = [sortedstat_p3 1]; 

% QOI_1: Maximum size of infectious population
h1=figure(1);
b = bar([sortedPRCC_1;SsortedPRCC_1]','FaceAlpha',.5);
% hold on
% plot(4,0,'xk','MarkerSize',24)
% hold off
legend('SEIR','SEIRS','fontsize',22)
legend boxoff
set(gca,'xticklabel',Ssortedparams,'fontsize',22)
ylabel('PRCC')
ylim([-1 1])

xtips = b(1).XEndPoints;
ytips = b(1).YEndPoints;

for j = 1:pars
    if (0.001 < sortedstat_p1(j)) && (sortedstat_p1(j) < 0.01)
        labels{j} = sprintf('**');
    elseif sortedstat_p1(j) < 0.001
        labels{j} = sprintf('*');
    else
        labels{j} = '';
    end
    if ytips(j) > 0
        text(xtips(j),ytips(j),labels(j),'fontsize',14,'HorizontalAlignment','center',...
            'VerticalAlignment','bottom')
    else
        text(xtips(j),ytips(j),labels(j),'fontsize',14,'HorizontalAlignment','center',...
            'VerticalAlignment','top')
    end
end


xtips = b(2).XEndPoints;
ytips = b(2).YEndPoints;
for j = 1:pars
    if (0.001 < Ssortedstat_p1(j)) && (Ssortedstat_p1(j) < 0.01)
        labels{j} = sprintf('**');
    elseif Ssortedstat_p1(j) < 0.001
        labels{j} = sprintf('*');
    else
        labels{j} = sprintf('%.2f', Ssortedstat_p1(j));
    end
    if ytips(j) > 0
        text(xtips(j),ytips(j),labels(j),'fontsize',14,'HorizontalAlignment','center',...
            'VerticalAlignment','bottom')
    else
        text(xtips(j),ytips(j),labels(j),'fontsize',14,'HorizontalAlignment','center',...
            'VerticalAlignment','top')
    end
end
saveas(h1,'QOI_1_SEIR_SEIRS.png')

% QOI_2: Infectious population size at 50 days
h2=figure(2);
b = bar([sortedPRCC_2;SsortedPRCC_2]','FaceAlpha',.5);
% hold on
% plot(4,0,'xk','MarkerSize',24)
% hold off
legend('SEIR','SEIRS','fontsize',22)
legend boxoff
set(gca,'xticklabel',Ssortedparams,'fontsize',22)
ylabel('PRCC')
ylim([-1 1])


xtips = b(1).XEndPoints;
ytips = b(1).YEndPoints;
for j = 1:pars
    if (0.001 < sortedstat_p2(j)) && (sortedstat_p2(j) < 0.01)
        labels{j} = sprintf('**');
    elseif sortedstat_p2(j) < 0.001
        labels{j} = sprintf('*');
    else
        labels{j} = '';
    end
    if ytips(j) > 0
        text(xtips(j),ytips(j),labels(j),'fontsize',14,'HorizontalAlignment','center',...
            'VerticalAlignment','bottom')
    else
        text(xtips(j),ytips(j),labels(j),'fontsize',14,'HorizontalAlignment','center',...
            'VerticalAlignment','top')
    end
end

xtips = b(2).XEndPoints;
ytips = b(2).YEndPoints;
for j = 1:pars
    if (0.001 < Ssortedstat_p2(j)) && (Ssortedstat_p2(j) < 0.01)
        labels{j} = sprintf('**');
    elseif Ssortedstat_p2(j) < 0.001
        labels{j} = sprintf('*');
    else
        labels{j} = sprintf('%.2f', Ssortedstat_p2(j));
    end
    if ytips(j) > 0
        text(xtips(j),ytips(j),labels(j),'fontsize',14,'HorizontalAlignment','center',...
            'VerticalAlignment','bottom')
    else
        text(xtips(j),ytips(j),labels(j),'fontsize',14,'HorizontalAlignment','center',...
            'VerticalAlignment','top')
    end
end
saveas(h2,'QOI_2_SEIR_SEIRS.png')

% QOI_3: Infectious population size at 2 years
h3=figure(3);
b = bar([sortedPRCC_3;SsortedPRCC_3]','FaceAlpha',.5);
set(gca,'xticklabel',Ssortedparams,'fontsize',22)
legend('SEIR','SEIRS','fontsize',22)
legend boxoff
ylabel('PRCC')
ylim([-1 1])


xtips = b(1).XEndPoints;
ytips = b(1).YEndPoints;
for j = 1:pars
    if (0.001 < sortedstat_p3(j)) && (sortedstat_p3(j) < 0.01)
        labels{j} = sprintf('**');
    elseif sortedstat_p3(j) < 0.001
        labels{j} = sprintf('*');
    else
        labels{j} = '';
    end
    if ytips(j) > 0
        text(xtips(j),ytips(j),labels(j),'fontsize',14,'HorizontalAlignment','center',...
            'VerticalAlignment','bottom')
    else
        text(xtips(j),ytips(j),labels(j),'fontsize',14,'HorizontalAlignment','center',...
            'VerticalAlignment','top')
    end
end

xtips = b(2).XEndPoints;
ytips = b(2).YEndPoints;
for j = 1:pars
    if (0.001 < Ssortedstat_p3(j)) && (Ssortedstat_p3(j) < 0.01)
        labels{j} = sprintf('**');
    elseif Ssortedstat_p3(j) < 0.001
        labels{j} = sprintf('*');
    else
        labels{j} = sprintf('%.2f', Ssortedstat_p3(j));
    end
    if ytips(j) > 0
        text(xtips(j),ytips(j),labels(j),'fontsize',14,'HorizontalAlignment','center',...
            'VerticalAlignment','bottom')
    else
        text(xtips(j),ytips(j),labels(j),'fontsize',14,'HorizontalAlignment','center',...
            'VerticalAlignment','top')
    end
end
saveas(h3,'QOI_3_SEIR_SEIRS.png')