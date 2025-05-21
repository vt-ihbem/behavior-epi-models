clear; close all; clc

%% LHS sampling
params = {'\alpha','\beta','\gamma','\tau_E','\tau_I','\tau_F','N','\tau_R'}; % parameter names
lb =     [   0        0.1      0        1         1       1     3         1]; % lower bound of parameters
ub =     [ 100          4      5        5        14     200     8       365]; % upper bound of parametersrt5f

samp = 10^4; % number of samples
pars = length(params); % number of parameters

load("LHS_scaled_sample.mat","LHS")

%% PRCC calculation

tmax=350*2;
IC = [1-(20/10^6), 10/10^6, 10/10^6, 0, 0];

% Adjust sample size here when troubleshooting; otherwise, set size1=samp
size1 = samp;
LHS = LHS(1:size1,1:8);

% Set zero matrices for gathering values for calculating quantities of interest (QOIs)
y50i=zeros(size1,1);
yendi=zeros(size1,1);
ymax=zeros(size1,1);

for j=1:size1 
    par.alpha = LHS(j,1);
    par.beta = LHS(j,2);
    par.gamma = LHS(j,3);
    par.tauE = LHS(j,4);
    par.tauI = LHS(j,5);
    par.tauF = LHS(j,6);
    par.N = LHS(j,7);
    par.tauR = LHS(j,8);

    options1=odeset('NonNegative',1:5,'RelTol',1e-9,'AbsTol',1e-12);
    [t,y]=ode45(@SEIRSb, 0:1:tmax, IC, options1, par);

    % Save max i value for PRCC
    ymax(j)=max(y(:,3));

    % Save all i values on day 50 for PRCC
    y50i(j)=y(50,3);

    % Save all i values after 2 years for PRCC
    yendi(j)=y(end,3);

end
hold off
%% Generate PRCCs and p-values for QOIs

QOI_1=ymax;
QOI_2=y50i;
QOI_3=yendi;

% Determine correlations
[rho1,p1] = partialcorr([LHS QOI_1],'type','Spearman'); 
[rho2,p2] = partialcorr([LHS QOI_2],'type','Spearman'); 
[rho3,p3] = partialcorr([LHS QOI_3],'type','Spearman');

% Correlations controlling for other parameters 
PRCC_1 = rho1(1:end-1,end)'; 
PRCC_2 = rho2(1:end-1,end)';
PRCC_3 = rho3(1:end-1,end)';

% Associated p-value
stat_p1 = p1(1:end-1,end)'; 
stat_p2 = p2(1:end-1,end)';
stat_p3 = p3(1:end-1,end)'; 

%% Plots: Check monotonicity
% Check monotonicity of relationships - do not include these plots in the final paper

% QOI_1
h1=figure(1);
for j = 1:pars
    subplot(2,4,j)
    semilogy(LHS(:,j),QOI_1,'.')
    xlabel(params{j})
    xlim([lb(j) ub(j)])
    set(gca,'fontsize',14)   
end
% Adjust axes in graph for parameter N
figure(1)
subplot(2,4,7)
xlim([10^lb(7) 10^ub(7)])
sgtitle('SEIRb: Max i - Check for monotonic relationships')
saveas(h1,'Mono_QOI_1_SEIRSb.png')

% QOI_2
h11=figure(11);
for j = 1:pars
    subplot(2,4,j)
    semilogy(LHS(:,j),QOI_2,'.')
    xlabel(params{j})
    xlim([lb(j) ub(j)])
    set(gca,'fontsize',14)   
end
% Adjust axes in graph for parameter N
figure(11)
subplot(2,4,7)
xlim([10^lb(7) 10^ub(7)])
sgtitle('SEIRb: i at 50 days - Check for monotonic relationships')
saveas(h11,'Mono_QOI_2_SEIRSb.png')

% QOI_3
h21=figure(21);
for j = 1:pars
    subplot(2,4,j)
    semilogy(LHS(:,j),QOI_3,'.')
    xlabel(params{j})
    xlim([lb(j) ub(j)])
    set(gca,'fontsize',14)   
end
% Adjust axes in graph for parameter N
figure(21)
subplot(2,4,7)
xlim([10^lb(7) 10^ub(7)]);
sgtitle('SEIRb: i at 2 years - Check for monotonic relationships')
saveas(h21,'Mono_QOI_3_SEIRSb.png')


%% Sort paramers here to match the order of those from the SEIRb model
load("LHS_sortOrder_SEIRb.mat","sortOrder")
sortOrder(8)=8;
sortedparams = params(sortOrder);

sortedPRCC_1 = PRCC_1(sortOrder);
sortedstat_p1 = stat_p1(sortOrder);


sortedPRCC_2 = PRCC_2(sortOrder);
sortedstat_p2 = stat_p2(sortOrder);


sortedPRCC_3 = PRCC_3(sortOrder);
sortedstat_p3 = stat_p3(sortOrder);

%% Plots: PRCC and p-values

% QOI_1
h2=figure(2);
b = bar(sortedPRCC_1,'FaceAlpha',.5);
set(gca,'xticklabel',sortedparams,'fontsize',22)
ylabel('PRCC')
ylim([-1 1])
%title('SEIRb: maximum size of infectious population') 

xtips = b.XEndPoints;
ytips = b.YEndPoints;
for j = 1:pars
    if (0.001 < sortedstat_p1(j)) && (sortedstat_p1(j) < 0.01)
        labels{j} = sprintf('*');
    elseif sortedstat_p1(j) < 0.001
        labels{j} = sprintf('**');
    else 
        labels{j} = sprintf('%.2f',sortedstat_p1(j));
    end
    if ytips(j) > 0
        text(xtips(j),ytips(j),labels(j),'fontsize',14,'HorizontalAlignment','center',...
            'VerticalAlignment','bottom')
    else
        text(xtips(j),ytips(j),labels(j),'fontsize',14,'HorizontalAlignment','center',...
            'VerticalAlignment','top')
    end
end
saveas(h2,'PRCC_QOI_1_SEIRSb.png')
saveas(h2,'PRCC_QOI_1_SEIRSb.fig')
saveas(h2,'PRCC_QOI_1_SEIRSb.pdf')


% QOI_2
h12=figure(12);
b = bar(sortedPRCC_2,'FaceAlpha',.5);
set(gca,'xticklabel',sortedparams,'fontsize',22)
ylabel('PRCC')
ylim([-1 1])
%title('SEIRb: infectious population size at 50 days')

xtips = b.XEndPoints;
ytips = b.YEndPoints;
for j = 1:pars
    if (0.001 < sortedstat_p2(j)) && (sortedstat_p2(j) < 0.01);
        labels{j} = sprintf('*');
    elseif sortedstat_p2(j) < 0.001
        labels{j} = sprintf('**');
    else 
        labels{j} = sprintf('%.2f',sortedstat_p2(j));
    end
    if ytips(j) > 0
        text(xtips(j),ytips(j),labels(j),'fontsize',14,'HorizontalAlignment','center',...
            'VerticalAlignment','bottom')
    else
        text(xtips(j),ytips(j),labels(j),'fontsize',14,'HorizontalAlignment','center',...
            'VerticalAlignment','top')
    end
end
saveas(h12,'PRCC_QOI_2_SEIRSb.png')
saveas(h12,'PRCC_QOI_2_SEIRSb.fig')
saveas(h12,'PRCC_QOI_2_SEIRSb.pdf')

% QOI_3
h22=figure(22);
b = bar(sortedPRCC_3,'FaceAlpha',.5);
set(gca,'xticklabel',sortedparams,'fontsize',22)
ylabel('PRCC')
ylim([-1 1])
%title('SEIRb: infectious population size at 2 years')

xtips = b.XEndPoints;
ytips = b.YEndPoints;
for j = 1:pars
    if (0.001 < sortedstat_p3(j)) &&  (sortedstat_p3(j) < 0.01)
        labels{j} = sprintf('*');
    elseif sortedstat_p3(j) < 0.001
        labels{j} = sprintf('**');
    else 
        labels{j} = sprintf('%.2f',sortedstat_p3(j));
    end
    if ytips(j) > 0
        text(xtips(j),ytips(j),labels(j),'fontsize',14,'HorizontalAlignment','center',...
            'VerticalAlignment','bottom')
    else
        text(xtips(j),ytips(j),labels(j),'fontsize',14,'HorizontalAlignment','center',...
            'VerticalAlignment','top')
    end
end
saveas(h22,'PRCC_QOI_3_SEIRSb.png')
saveas(h22,'PRCC_QOI_3_SEIRSb.fig')
saveas(h22,'PRCC_QOI_3_SEIRSb.pdf')

save("LHS_PRCC_SEIRSb_output.mat",...
     'sortedPRCC_1','sortedPRCC_2','sortedPRCC_3','sortedstat_p1',...
     'sortedstat_p2','sortedstat_p3','sortedparams','pars')


%{ %% Quantity of Interest 1 (QOI_1): max i (infectious population size)
% QOI_1=ymax;
% 
% [rho1,p1] = partialcorr([LHS QOI_1],'type','Spearman'); % determine correlations
% PRCC_1 = rho1(1:end-1,end)'; % correlations controling for other parameters 
% stat_p1 = p1(1:end-1,end)'; % associated p-value
% 
% 
% % Check monotonicity of relationships - do not include this plot in the final paper
% 
% h1=figure(1);
% for j = 1:pars-1
%     subplot(2,4,j)
%     semilogy(LHS(:,j),QOI_1,'.')
%     xlabel(params{j})
%     xlim([lb(j) ub(j)])
%     set(gca,'fontsize',14)   
% end
% % Adjust axes in graph for parameter N
% figure(1)
% subplot(2,4,8)
% xlim([10^lb(8) 10^ub(8)])
% sgtitle('SEIRSb: Max i - Check for monotonic relationships')
% saveas(h1,'Mono_QOI_1_SEIRSb.png')
% 
% 
% % Sort parameters using order from SEIRb LHS/PRCC analysis
% %sortOrder = [5 2 6 1 3 8 4 7];
% load("LHS_sortOrder_SEIRb.mat","sortOrder");
% sortedPRCC_1 = PRCC_1(sortOrder);
% 
% % Sort paramers here to order from most to least significant impact on max i; this order will remain the same for the other QOIs
% sortedparams = params(sortOrder);
% sortedstat_p1 = stat_p1(sortOrder);
% 
% h2=figure(2);
% b = bar(sortedPRCC_1,'FaceAlpha',.5);
% set(gca,'xticklabel',sortedparams,'fontsize',22)
% ylabel('PRCC')
% ylim([-1 1])
% %title('SEIRSb: maximum size of infectious population')
% 
% xtips = b.XEndPoints;
% ytips = b.YEndPoints;
% for j = 1:pars
%     if 0.001 < sortedstat_p1(j) < 0.01
%         labels{j} = sprintf('*',sortedstat_p1(j));
%     elseif sortedstat_p1(j) < 0.001
%         labels{j} = sprintf('**',sortedstat_p1(j));
%     else 
%         labels{j} = sprintf('%.2f',sortedstat_p1(j));
%     end
%     if ytips(j) > 0
%         text(xtips(j),ytips(j),labels(j),'fontsize',14,'HorizontalAlignment','center',...
%             'VerticalAlignment','bottom')
%     else
%         text(xtips(j),ytips(j),labels(j),'fontsize',14,'HorizontalAlignment','center',...
%             'VerticalAlignment','top')
%     end
% end
% saveas(h2,'PRCC_QOI_1_SEIRSb.png')
% saveas(h2,'PRCC_QOI_1_SEIRSb.fig')
% 
% 
% %% Quantity of Interest 2 (QOI_2): i at 50 days
% 
% QOI_2=y50i;
% 
% [rho2,p2] = partialcorr([LHS QOI_2],'type','Spearman'); % determine correlations
% 
% PRCC_2 = rho2(1:end-1,end)'; % correlations controling for other parameters 
% stat_p2 = p2(1:end-1,end)'; % associated p-value
% 
% % Check monotonicity of relationships - do not include this plot in the final paper
% h11=figure(11);
% for j = 1:pars
%     subplot(2,4,j)
%     semilogy(LHS(:,j),QOI_2,'.')
%     xlabel(params{j})
%     xlim([lb(j) ub(j)])
%     set(gca,'fontsize',14)   
% end
% % Adjust axes in graph for parameter N
% figure(11)
% subplot(2,4,8)
% xlim([10^lb(8) 10^ub(8)])
% sgtitle('SEIRSb: i at 50 days - Check for monotonic relationships')
% saveas(h11,'Mono_QOI_2_SEIRSb.png')
% 
% sortedPRCC_2 = PRCC_2(sortOrder);
% sortedstat_p2 = stat_p2(sortOrder);
% 
% h12=figure(12);
% b = bar(sortedPRCC_2,'FaceAlpha',.5);
% set(gca,'xticklabel',sortedparams,'fontsize',22)
% ylabel('PRCC')
% ylim([-1 1])
% %title('SEIRSb: infectious population size at 50 days')
% 
% xtips = b.XEndPoints;
% ytips = b.YEndPoints;
% for j = 1:pars
%     if 0.001 < sortedstat_p2(j) < 0.01
%         labels{j} = sprintf('*',sortedstat_p2(j));
%     elseif sortedstat_p2(j) < 0.001
%         labels{j} = sprintf('**',sortedstat_p2(j));
%     else 
%         labels{j} = sprintf('%.2f',sortedstat_p2(j));
%     end
%     if ytips(j) > 0
%         text(xtips(j),ytips(j),labels(j),'fontsize',14,'HorizontalAlignment','center',...
%             'VerticalAlignment','bottom')
%     else
%         text(xtips(j),ytips(j),labels(j),'fontsize',14,'HorizontalAlignment','center',...
%             'VerticalAlignment','top')
%     end
% end
% saveas(h12,'PRCC_QOI_2_SEIRSb.png')
% saveas(h12,'PRCC_QOI_2_SEIRSb.fig')
% 
% 
% %% Quantity of Interest 3 (QOI_3): i at 2 years
% 
% QOI_3=yendi;
% 
% % For QOI_3:
% [rho3,p3] = partialcorr([LHS QOI_3],'type','Spearman'); % determine correlations
% 
% PRCC_3 = rho3(1:end-1,end)'; % correlations controling for other parameters 
% stat_p3 = p3(1:end-1,end)'; % associated p-value
% 
% % Check monotonicity of relationships - do not include this plot in the final paper
% h21 = figure(21);
% for j = 1:pars
%     subplot(2,4,j)
%     semilogy(LHS(:,j),QOI_3,'.')
%     xlabel(params{j})
%     xlim([lb(j) ub(j)])
%     set(gca,'fontsize',14)   
% end
% % Adjust axes in graph for parameter N
% figure(21)
% subplot(2,4,8)
% xlim([10^lb(8) 10^ub(8)])
% sgtitle('SEIRSb: i at 2 years - Check for monotonic relationships')
% saveas(h21,'Mono_QOI_3_SEIRSb.png')
% 
% sortedPRCC_3 = PRCC_3(sortOrder);
% sortedstat_p3 = stat_p3(sortOrder);
% 
% h22=figure(22);
% b = bar(sortedPRCC_3,'FaceAlpha',.5);
% set(gca,'xticklabel',sortedparams,'fontsize',22)
% ylabel('PRCC')
% ylim([-1 1])
% %title('SEIRSb: infectious population size at 2 years')
% 
% xtips = b.XEndPoints;
% ytips = b.YEndPoints;
% for j = 1:pars
%     if 0.001 < sortedstat_p3(j) < 0.01
%         labels{j} = sprintf('*',sortedstat_p3(j));
%     elseif sortedstat_p3(j) < 0.001
%         labels{j} = sprintf('**',sortedstat_p3(j));
%     else 
%         labels{j} = sprintf('%.2f',sortedstat_p3(j));
%     end
%     if ytips(j) > 0
%         text(xtips(j),ytips(j),labels(j),'fontsize',14,'HorizontalAlignment','center',...
%             'VerticalAlignment','bottom')
%     else
%         text(xtips(j),ytips(j),labels(j),'fontsize',14,'HorizontalAlignment','center',...
%             'VerticalAlignment','top')
%     end
% end
% saveas(h22,'PRCC_QOI_3_SEIRSb.png')
%} saveas(h22,'PRCC_QOI_3_SEIRSb.fig')

function dy = SEIRSb(~,y,p)
s=y(1); e=y(2); i=y(3); r=y(4); f=y(5);
dy=zeros(5,1);

dy(1)=-p.beta/((1+p.alpha*p.N*f)^p.gamma)*s*i+r/p.tauR;
dy(2)= p.beta/((1+p.alpha*p.N*f)^p.gamma)*s*i-e/p.tauE;
dy(3)= e/p.tauE-i/p.tauI;
dy(4)= i/p.tauI-r/p.tauR;
dy(5)= (i-f)/p.tauF;
end
