clearvars; clc; close all; 

p.beta  = 0.7;
p.tauE  = 5;
p.tauI  = 10;
p.tauR  = 90;

% Perturb beta down
p1b.beta = 0.1; 
p1b.tauE=p.tauE;
p1b.tauI=p.tauI;
p1b.tauR = p.tauR;

p2b.beta = 0.4;
p2b.tauE=p.tauE;
p2b.tauI=p.tauI;
p2b.tauR = p.tauR;

% Perturb beta up
p3b.beta = 1;
p3b.tauE=p.tauE;
p3b.tauI=p.tauI;
p3b.tauR = p.tauR;

p4b.beta = 1.3;
p4b.tauE=p.tauE;
p4b.tauI=p.tauI;
p4b.tauR = p.tauR;

% Perturb tauR down
p1tauR.beta = p.beta; 
p1tauR.tauE=p.tauE;
p1tauR.tauI=p.tauI;
p1tauR.tauR = 30;

p2tauR.beta = p.beta; 
p2tauR.tauE=p.tauE;
p2tauR.tauI=p.tauI;
p2tauR.tauR = 60;

% Perturb tauR up
p3tauR.beta = p.beta; 
p3tauR.tauE=p.tauE;
p3tauR.tauI=p.tauI;
p3tauR.tauR = 120;

p4tauR.beta = p.beta; 
p4tauR.tauE=p.tauE;
p4tauR.tauI=p.tauI;
p4tauR.tauR = 150;


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

% Solution w/tauR perturbation
[t1tauR,y1tauR]=ode45(@sens, [0:0.01:tmax], IC, [], p1tauR);
[t2tauR,y2tauR]=ode45(@sens, [0:0.01:tmax], IC, [], p2tauR);
[t3tauR,y3tauR]=ode45(@sens, [0:0.01:tmax], IC, [], p3tauR);
[t4tauR,y4tauR]=ode45(@sens, [0:0.01:tmax], IC, [], p4tauR);


% Plot for beta
h1=figure(1);
hold on
plot(t1b,y1b(:,3),'r','LineWidth',2)
plot(t2b,y2b(:,3),'b','LineWidth',2)
plot(t,y(:,3),'--g','LineWidth',2)
plot(t3b,y3b(:,3),'k','LineWidth',2)
plot(t4b,y4b(:,3),'m','LineWidth',2)
xlim([0,tmax]);
xlabel('Time (days)','FontSize',24)
ylabel('Infectious population size (i)','FontSize',24)
legend(sprintf('\\beta = %g',p1b.beta),sprintf('\\beta = %g',p2b.beta),...
        sprintf('\\beta = %g',p.beta),sprintf('\\beta = %g',p3b.beta),sprintf('\\beta = %g',p4b.beta),...
        'FontSize',22)
legend boxoff
ax = gca;
ax.XAxis.FontSize = 24;
ax.YAxis.FontSize = 24;
saveas(h1,'Change_beta_SEIRS.png');

% Plot for tauR
h2=figure(2);
hold on
plot(t1tauR,y1tauR(:,3),'r','LineWidth',2)
plot(t2tauR,y2tauR(:,3),'b','LineWidth',2)
plot(t,y(:,3),'--g','LineWidth',2)
plot(t3tauR,y3tauR(:,3),'k','LineWidth',2)
plot(t4tauR,y4tauR(:,3),'m','LineWidth',2)
xlim([0,tmax]);
ylim([0,0.5]);
xlabel('Time (days)','FontSize',24)
ylabel('Infectious population size (i)','FontSize',24)
legend(sprintf('\\tau_R = %g',p1tauR.tauR),sprintf('\\tau_R = %g',p2tauR.tauR),sprintf('\\tau_R = %g',p.tauR),...
       sprintf('\\tau_R = %g',p3tauR.tauR),sprintf('\\tau_R = %g',p4tauR.tauR),'FontSize',22)
legend boxoff
ax = gca;
ax.XAxis.FontSize = 24;
ax.YAxis.FontSize = 24;
saveas(h2,'Change_tauR_SEIRS.png');


function dy = sens(~,y,p)
S=y(1); E=y(2); I=y(3); R=y(4);
dy = zeros(4,1);

% Differential equations for S, E, I, R
dy(1) = -p.beta*S*I+R/p.tauR;
dy(2) =  p.beta*S*I-E/p.tauE;
dy(3) =  E/p.tauE-I/p.tauI;
dy(4) =  I/p.tauI-R/p.tauR;

end



