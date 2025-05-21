clear; close all; clc

%% LHS sampling
params = {'\beta','\tau_E','\tau_I','\tau_R'}; % parameter names
lb =     [0.1       1         1            1]; % lower bound of parameters
ub =     [4         5        14          365]; % upper bound of parameters

samp = 10^4; % number of samples
pars = length(params);

load("LHS_raw.mat","LHS_raw");

% Pull relevant params for SEIR model
LHS_raw1(:,1) = LHS_raw(:,2); 
LHS_raw1(:,2) = LHS_raw(:,4); 
LHS_raw1(:,3) = LHS_raw(:,5); 
LHS_raw1(:,4) = LHS_raw(:,8);

% Scale parameters to range of interest
LHS = LHS_raw1.*(ub-lb)+lb; 
save("LHS_scaled_sample.mat","LHS") % save LHS output to allow for exact reproduction of results

%% PRCC calculation

tmax=350*2;
IC = [1-(20/10^6), 10/10^6, 10/10^6, 0];

% Adjust sample size here when troubleshooting; otherwise, set size1=samp
size1 = samp;
LHS = LHS(1:size1,1:4);

% Set zero matrices for gathering values for calculating quantities of interest (QOIs)
y50i=zeros(size1,1);
yendi=zeros(size1,1);
ymax=zeros(size1,1);

%% Solver

for j=1:size1 
    par.beta = LHS(j,1);
    par.tauE = LHS(j,2);
    par.tauI = LHS(j,3);
    par.tauR = LHS(j,4);

    options1=odeset('NonNegative',1:4,'RelTol',1e-9,'AbsTol',1e-12);
    [t,y]=ode45(@SEIRS, 0:1:tmax, IC, options1, par);

    % Save max i value for PRCC
    ymax(j)=max(y(:,3));

    % Save all i values on day 50 for PRCC
    y50i(j)=y(50,3);

    % Save all i values after 2 years for PRCC
    yendi(j)=y(end,3);

end

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
    subplot(1,4,j)
    plot(LHS(:,j),QOI_1,'.')
    xlabel(params{j})
    xlim([lb(j) ub(j)])
    set(gca,'fontsize',14)   
end
sgtitle('SEIRS: Max i - Check for monotonic relationships')
saveas(h1,'Mono_QOI_1_SEIRS.png')

% QOI_2
h11=figure(11);
for j = 1:pars
    subplot(1,4,j)
    plot(LHS(:,j),QOI_2,'.')
    xlabel(params{j})
    xlim([lb(j) ub(j)])
    set(gca,'fontsize',14)   
end
sgtitle('SEIRS: i at 50 days - Check for monotonic relationships')
saveas(h11,'Mono_QOI_2_SEIRS.png')

% QOI_3
h21=figure(21);
for j = 1:pars
    subplot(1,4,j)
    plot(LHS(:,j),QOI_3,'.')
    xlabel(params{j})
    xlim([lb(j) ub(j)])
    set(gca,'fontsize',14)   
end
sgtitle('SEIRS: i at 2 years - Check for monotonic relationships')
saveas(h21,'Mono_QOI_3_SEIRS.png')

%% Sort for PRCC plotting

% Sort parameters here to order from most to least significant impact on max i;
% this order will remain the same for the other QOIs

% Sort elements in order of SEIRb: tauI, beta, tauE, tauR
K=[2 3 1 4];
[~,sortOrder] = sort(K);
sortedparams = params(sortOrder)

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
%title('SEIRS: maximum size of infectious population') 

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
saveas(h2,'PRCC_QOI_1_SEIRS.png')
saveas(h2,'PRCC_QOI_1_SEIRS.fig')
saveas(h2,'PRCC_QOI_1_SEIRS.pdf')


% QOI_2
h12=figure(12);
b = bar(sortedPRCC_2,'FaceAlpha',.5);
set(gca,'xticklabel',sortedparams,'fontsize',22)
ylabel('PRCC')
ylim([-1 1])
%title('SEIRS: infectious population size at 50 days')

xtips = b.XEndPoints;
ytips = b.YEndPoints;
for j = 1:pars
    if (0.001 < sortedstat_p2(j)) && (sortedstat_p2(j) < 0.01)
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
saveas(h12,'PRCC_QOI_2_SEIRS.png')
saveas(h12,'PRCC_QOI_2_SEIRS.fig')
saveas(h12,'PRCC_QOI_2_SEIRS.pdf')

% QOI_3
h22=figure(22);
b = bar(sortedPRCC_3,'FaceAlpha',.5);
set(gca,'xticklabel',sortedparams,'fontsize',22)
ylabel('PRCC')
ylim([-1 1])
%title('SEIRS: infectious population size at 2 years')

xtips = b.XEndPoints;
ytips = b.YEndPoints;
for j = 1:pars
    if (0.001 < sortedstat_p3(j)) &&  (sortedstat_p3(j) < 0.01)
        labels{j} = sprintf('*');
    elseif sortedstat_p3(j) < 0.001
        labels{j} = sprintf('**');
    else 
        labels{j} = sprintf('%.2f', sortedstat_p3(j));
    end
    if ytips(j) > 0
        text(xtips(j),ytips(j),labels(j),'fontsize',14,'HorizontalAlignment','center',...
            'VerticalAlignment','bottom')
    else
        text(xtips(j),ytips(j),labels(j),'fontsize',14,'HorizontalAlignment','center',...
            'VerticalAlignment','top')
    end
end
saveas(h22,'PRCC_QOI_3_SEIRS.png')
saveas(h22,'PRCC_QOI_3_SEIRS.fig')
saveas(h22,'PRCC_QOI_3_SEIRS.pdf')

save("LHS_PRCC_SEIRS_output.mat",...
     'sortedPRCC_1','sortedPRCC_2','sortedPRCC_3','sortedstat_p1',...
     'sortedstat_p2','sortedstat_p3','sortedparams','pars')

%% Function

function dy = SEIRS(~,y,p)
s=y(1); e=y(2); i=y(3); r=y(4);
dy=zeros(4,1);

dy(1)=-p.beta*s*i+r/p.tauR;
dy(2)= p.beta*s*i-e/p.tauE;
dy(3)= e/p.tauE-i/p.tauI;
dy(4)= i/p.tauI-r/p.tauR;

end
