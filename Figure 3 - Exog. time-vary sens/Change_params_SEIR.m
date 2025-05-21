clearvars; close all; clc

p.beta  = 0.7;
p.tauE  = 5;
p.tauI  = 10;

% Perturb beta down
p1b.beta = 0.1; 
p1b.tauE=p.tauE;
p1b.tauI=p.tauI;

p2b.beta = 0.4;
p2b.tauE=p.tauE;
p2b.tauI=p.tauI;

% Perturb beta up
p3b.beta = 1;
p3b.tauE=p.tauE;
p3b.tauI=p.tauI;

p4b.beta = 1.3;
p4b.tauE=p.tauE;
p4b.tauI=p.tauI;

S0 = 1-(20/10^6);
E0 = 10/10^6;
I0 = 10/10^6;
R0 = 0;

IC = [S0, E0, I0, R0];
tmax=350;

% Solution w/no parameter perturbation
[t,y] = ode15s(@sens, [0:0.01:tmax], IC, [], p);

% Solution w/beta perturbation
[t1b,y1b]=ode45(@sens, [0:0.01:tmax], IC, [], p1b);
[t2b,y2b]=ode45(@sens, [0:0.01:tmax], IC, [], p2b);
[t3b,y3b]=ode45(@sens, [0:0.01:tmax], IC, [], p3b);
[t4b,y4b]=ode45(@sens, [0:0.01:tmax], IC, [], p4b);
% m1 = max(y1b(:,3));
% m2 = max(y2b(:,3)); 
% m3 = max(y3b(:,3));
% m4 = max(y4b(:,3));
%m = max([max(y1b(:,3)) max(y2b(:,3)) max(y3b(:,3)) max(y4b(:,4))])
%m = max([m1, m2, m3, m4]);

% Plot for beta
h1=figure(1);
hold on
plot(t1b,y1b(:,3),'r','LineWidth',2)
plot(t2b,y2b(:,3),'b','LineWidth',2)
plot(t,y(:,3),'--g','LineWidth',2)
plot(t3b,y3b(:,3),'k','LineWidth',2)
plot(t4b,y4b(:,3),'m','LineWidth',2)
xlim([0,tmax]);
ylim([0,0.5]);
xlabel('Time (days)','FontSize',24)
ylabel('Infectious population size (i)','FontSize',24)
legend(sprintf('\\beta = %g',p1b.beta),sprintf('\\beta = %g',p2b.beta),...
        sprintf('\\beta = %g',p.beta),sprintf('\\beta = %g',p3b.beta),sprintf('\\beta = %g',p4b.beta),...
        'FontSize',22)
legend boxoff
ax = gca;
ax.XAxis.FontSize = 24;
ax.YAxis.FontSize = 24;
saveas(h1,'Change_beta_SEIR.png');


function dy = sens(~,y,p)
S=y(1); E=y(2); I=y(3); R=y(4);
dy = zeros(4,1);

% Differential equations for S, E, I, R
dy(1) = -p.beta*S*I;
dy(2) =  p.beta*S*I-E/p.tauE;
dy(3) =  E/p.tauE-I/p.tauI;
dy(4) =  I/p.tauI;

end


