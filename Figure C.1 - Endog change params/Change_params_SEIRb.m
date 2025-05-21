clearvars; close all; clc

%% Original parameters
p.beta  = 0.7;
p.a = 100/3000*10^6;
p.gamma = 2;
p.tauE  = 5;
p.tauI  = 10;
p.tauF  = 20;

%% Change beta
% Perturb beta down
p1b.beta = 0.1; 
p1b.a = p.a;
p1b.gamma=p.gamma;
p1b.tauE=p.tauE;
p1b.tauI=p.tauI;
p1b.tauF = p.tauF;

p2b.beta = 0.4;
p2b.a = p.a;
p2b.gamma=p.gamma;
p2b.tauE=p.tauE;
p2b.tauI=p.tauI;
p2b.tauF = p.tauF;

% Perturb beta up
p3b.beta = 1;
p3b.a = p.a;
p3b.gamma=p.gamma;
p3b.tauE=p.tauE;
p3b.tauI=p.tauI;
p3b.tauF = p.tauF;

p4b.beta = 1.3;
p4b.a = p.a;
p4b.gamma=p.gamma;
p4b.tauE=p.tauE;
p4b.tauI=p.tauI;
p4b.tauF = p.tauF;


%% Change gamma
% Perturb gamma down
p1g.beta = p.beta; 
p1g.a = p.a;
p1g.gamma= 1;
p1g.tauE=p.tauE;
p1g.tauI=p.tauI;
p1g.tauF = p.tauF;

% Perturb gamma up
p2g.beta = p.beta; 
p2g.a = p.a;
p2g.gamma= 3;
p2g.tauE=p.tauE;
p2g.tauI=p.tauI;
p2g.tauF = p.tauF;

p3g.beta = p.beta; 
p3g.a = p.a;
p3g.gamma = 5;
p3g.tauE = p.tauE;
p3g.tauI = p.tauI;
p3g.tauF = p.tauF;

p4g.beta = p.beta; 
p4g.a = p.a;
p4g.gamma= 9;
p4g.tauE=p.tauE;
p4g.tauI=p.tauI;
p4g.tauF = p.tauF;


%% Change a
% Perturb a down
p1a.beta = p.beta; 
p1a.a = p.a -2*(.01*10^6);
p1a.gamma=p.gamma;
p1a.tauE=p.tauE;
p1a.tauI=p.tauI;
p1a.tauF = p.tauF;

p2a.beta = p.beta; 
p2a.a = p.a-.01*10^6;
p2a.gamma=p.gamma;
p2a.tauE=p.tauE;
p2a.tauI=p.tauI;
p2a.tauF = p.tauF;

% Perturb a up
p3a.beta = p.beta; 
p3a.a = p.a+.01*10^6;
p3a.gamma=p.gamma;
p3a.tauE=p.tauE;
p3a.tauI=p.tauI;
p3a.tauF = p.tauF;

p4a.beta = p.beta; 
p4a.a = p.a+2*(.01*10^6);
p4a.gamma=p.gamma;
p4a.tauE=p.tauE;
p4a.tauI=p.tauI;
p4a.tauF = p.tauF;


%% Change tauE
% Perturb tauE down
p1tauE.beta = p.beta; 
p1tauE.a = p.a;
p1tauE.gamma=p.gamma;
p1tauE.tauE=p.tauE-4;
p1tauE.tauI=p.tauI;
p1tauE.tauF = p.tauF;

p2tauE.beta = p.beta; 
p2tauE.a = p.a;
p2tauE.gamma=p.gamma;
p2tauE.tauE=p.tauE-3;
p2tauE.tauI=p.tauI;
p2tauE.tauF = p.tauF;

p3tauE.beta = p.beta; 
p3tauE.a = p.a;
p3tauE.gamma=p.gamma;
p3tauE.tauE=p.tauE-2;
p3tauE.tauI=p.tauI;
p3tauE.tauF = p.tauF;

p4tauE.beta = p.beta; 
p4tauE.a = p.a;
p4tauE.gamma=p.gamma;
p4tauE.tauE=p.tauE-1;
p4tauE.tauI=p.tauI;
p4tauE.tauF = p.tauF;

%% Change tauI
% Perturb tauI down
p1tauI.beta = p.beta; 
p1tauI.a = p.a;
p1tauI.gamma=p.gamma;
p1tauI.tauE=p.tauE;
p1tauI.tauI=p.tauI-4;
p1tauI.tauF = p.tauF;

p2tauI.beta = p.beta; 
p2tauI.a = p.a;
p2tauI.gamma=p.gamma;
p2tauI.tauE=p.tauE;
p2tauI.tauI=p.tauI-2;
p2tauI.tauF = p.tauF;

% Perturb tauI up
p3tauI.beta = p.beta; 
p3tauI.a = p.a;
p3tauI.gamma=p.gamma;
p3tauI.tauE=p.tauE;
p3tauI.tauI=p.tauI+2;
p3tauI.tauF = p.tauF;

p4tauI.beta = p.beta; 
p4tauI.a = p.a;
p4tauI.gamma=p.gamma;
p4tauI.tauE=p.tauE;
p4tauI.tauI=p.tauI+4;
p4tauI.tauF = p.tauF;


%% Change tauF
% Perturb tauF down
p1tauF.beta = p.beta; 
p1tauF.a = p.a;
p1tauF.gamma=p.gamma;
p1tauF.tauE=p.tauE;
p1tauF.tauI=p.tauI;
p1tauF.tauF = p.tauF-10;

% Perturb tauF up
p2tauF.beta = p.beta; 
p2tauF.a = p.a;
p2tauF.gamma=p.gamma;
p2tauF.tauE=p.tauE;
p2tauF.tauI=p.tauI;
p2tauF.tauF = p.tauF+10;

p3tauF.beta = p.beta; 
p3tauF.a = p.a;
p3tauF.gamma=p.gamma;
p3tauF.tauE=p.tauE;
p3tauF.tauI=p.tauI;
p3tauF.tauF = p.tauF+20;

p4tauF.beta = p.beta; 
p4tauF.a = p.a;
p4tauF.gamma=p.gamma;
p4tauF.tauE=p.tauE;
p4tauF.tauI=p.tauI;
p4tauF.tauF = p.tauF+30;

%% Initial conditions and solvers
S0 = 1-(20/10^6);
E0 = 10/10^6;
I0 = 10/10^6;
R0 = 0;
F0 = 0;

IC = [S0, E0, I0, R0, F0];
tmax=350;

% Solution w/no parameter perturbation
[t,y] = ode15s(@sens, [0:0.01:tmax], IC, [], p);

% Solution w/beta perturbation
[t1b,y1b]=ode45(@sens, [0:0.01:tmax], IC, [], p1b);
[t2b,y2b]=ode45(@sens, [0:0.01:tmax], IC, [], p2b);
[t3b,y3b]=ode45(@sens, [0:0.01:tmax], IC, [], p3b);
[t4b,y4b]=ode45(@sens, [0:0.01:tmax], IC, [], p4b);

% Solution w/gamma perturbation
[t1g,y1g]=ode45(@sens, [0:0.01:tmax], IC, [], p1g);
[t2g,y2g]=ode45(@sens, [0:0.01:tmax], IC, [], p2g);
[t3g,y3g]=ode45(@sens, [0:0.01:tmax], IC, [], p3g);
[t4g,y4g]=ode45(@sens, [0:0.01:tmax], IC, [], p4g);

% Solution w/a perturbation
[t1a,y1a]=ode45(@sens, [0:0.01:tmax], IC, [], p1a);
[t2a,y2a]=ode45(@sens, [0:0.01:tmax], IC, [], p2a);
[t3a,y3a]=ode45(@sens, [0:0.01:tmax], IC, [], p3a);
[t4a,y4a]=ode45(@sens, [0:0.01:tmax], IC, [], p4a);

% Solution w/tauE perturbation
[t1tauE,y1tauE]=ode45(@sens, [0:0.01:tmax], IC, [], p1tauE);
[t2tauE,y2tauE]=ode45(@sens, [0:0.01:tmax], IC, [], p2tauE);
[t3tauE,y3tauE]=ode45(@sens, [0:0.01:tmax], IC, [], p3tauE);
[t4tauE,y4tauE]=ode45(@sens, [0:0.01:tmax], IC, [], p4tauE);

% Solution w/tauI perturbation
[t1tauI,y1tauI]=ode45(@sens, [0:0.01:tmax], IC, [], p1tauI);
[t2tauI,y2tauI]=ode45(@sens, [0:0.01:tmax], IC, [], p2tauI);
[t3tauI,y3tauI]=ode45(@sens, [0:0.01:tmax], IC, [], p3tauI);
[t4tauI,y4tauI]=ode45(@sens, [0:0.01:tmax], IC, [], p4tauI);

% Solution w/tauF perturbation
[t1tauF,y1tauF]=ode45(@sens, [0:0.01:tmax], IC, [], p1tauF);
[t2tauF,y2tauF]=ode45(@sens, [0:0.01:tmax], IC, [], p2tauF);
[t3tauF,y3tauF]=ode45(@sens, [0:0.01:tmax], IC, [], p3tauF);
[t4tauF,y4tauF]=ode45(@sens, [0:0.01:tmax], IC, [], p4tauF);


%% Plots
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
        sprintf('\\beta = %g',p.beta),sprintf('\\beta = %g',p3b.beta),...
        sprintf('\\beta = %g',p4b.beta),'FontSize',24)
legend boxoff
ax = gca;
ax.XAxis.FontSize = 24;
ax.YAxis.FontSize = 24;
saveas(h1,'Change_beta_SEIRb.png');

% Plot for gamma
h2=figure(2);
hold on
plot(t1g,y1g(:,3),'r','LineWidth',2)
plot(t,y(:,3),'--g','LineWidth',2)
plot(t2g,y2g(:,3),'b','LineWidth',2)
plot(t3g,y3g(:,3),'k','LineWidth',2)
plot(t4g,y4g(:,3),'m','LineWidth',2)
xlim([0,tmax]);
xlabel('Time (days)','FontSize',24)
ylabel('Infectious population size (i)','FontSize',24)
legend(sprintf('\\gamma = %g',p1g.gamma),sprintf('\\gamma = %g',p.gamma),...
        sprintf('\\gamma = %g',p2g.gamma),sprintf('\\gamma = %g',p3g.gamma),...
        sprintf('\\gamma = %g',p4g.gamma),'FontSize',24)
legend boxoff
ax = gca;
ax.XAxis.FontSize = 24;
ax.YAxis.FontSize = 24;
saveas(h2,'Change_gamma_SEIRb.png');

% Plot for a
h3=figure(3);
hold on
plot(t1a,y1a(:,3),'r','LineWidth',2)
plot(t2a,y2a(:,3),'b','LineWidth',2)
plot(t,y(:,3),'--g','LineWidth',2)
plot(t3a,y3a(:,3),'k','LineWidth',2)
plot(t4a,y4a(:,3),'m','LineWidth',2)
xlim([0,tmax]);
xlabel('Time (days)','FontSize',24)
ylabel('Infectious population size (i)','FontSize',24)
legend(sprintf('a = %g',p1a.a),sprintf('a = %g',p2a.a),sprintf('a = %g',p.a),...
        sprintf('a = %g',p3a.a),sprintf('a = %g',p4a.a),'FontSize',24)
legend boxoff
ax = gca;
ax.XAxis.FontSize = 24;
ax.YAxis.FontSize = 24;
saveas(h3,'Change_a_SEIRb.png');

% Plot for tauE
h4=figure(4);
hold on
plot(t1tauE,y1tauE(:,3),'r','LineWidth',2)
plot(t2tauE,y2tauE(:,3),'b','LineWidth',2)
plot(t3tauE,y3tauE(:,3),'k','LineWidth',2)
plot(t4tauE,y4tauE(:,3),'m','LineWidth',2)
plot(t,y(:,3),'--g','LineWidth',2)
xlim([0,tmax]);
xlabel('Time (days)','FontSize',24)
ylabel('Infectious population size (i)','FontSize',24)
legend(sprintf('\\tau_E = %g',p1tauE.tauE),sprintf('\\tau_E = %g',p2tauE.tauE),...
        sprintf('\\tau_E = %g',p3tauE.tauE),sprintf('\\tau_E = %g',p4tauE.tauE),sprintf('\\tau_E = %g',p.tauE),'FontSize',24)
legend boxoff
ax = gca;
ax.XAxis.FontSize = 24;
ax.YAxis.FontSize = 24;
saveas(h4,'Change_tauE_SEIRb.png');

% Plot for tauI
h5=figure(5);
hold on
plot(t1tauI,y1tauI(:,3),'r','LineWidth',2)
plot(t2tauI,y2tauI(:,3),'b','LineWidth',2)
plot(t,y(:,3),'--g','LineWidth',2)
plot(t3tauI,y3tauI(:,3),'k','LineWidth',2)
plot(t4tauI,y4tauI(:,3),'m','LineWidth',2)
xlim([0,tmax]);
xlabel('Time (days)','FontSize',24)
ylabel('Infectious population size (i)','FontSize',24)
legend(sprintf('\\tau_I = %g',p1tauI.tauI),sprintf('\\tau_I = %g',p2tauI.tauI),sprintf('\\tau_I = %g',p.tauI),...
        sprintf('\\tau_I = %g',p3tauI.tauI),sprintf('\\tau_I = %g',p4tauI.tauI),'FontSize',24)
legend boxoff
ax = gca;
ax.XAxis.FontSize = 24;
ax.YAxis.FontSize = 24;
saveas(h5,'Change_tauI_SEIRb.png');

% Plot for tauF
h6=figure(6);
hold on
plot(t1tauF,y1tauF(:,3),'r','LineWidth',2)
plot(t,y(:,3),'--g','LineWidth',2)
plot(t2tauF,y2tauF(:,3),'b','LineWidth',2)
plot(t3tauF,y3tauF(:,3),'k','LineWidth',2)
plot(t4tauF,y4tauF(:,3),'m','LineWidth',2)
xlim([0,tmax]);
xlabel('Time (days)','FontSize',24)
ylabel('Infectious population size (i)','FontSize',24)
legend(sprintf('\\tau_F = %g',p1tauF.tauF),sprintf('\\tau_F = %g',p.tauF),sprintf('\\tau_F = %g',p2tauF.tauF),...
        sprintf('\\tau_F = %g',p3tauF.tauF),sprintf('\\tau_F = %g',p4tauF.tauF),'FontSize',24)
legend boxoff
ax = gca;
ax.XAxis.FontSize = 24;
ax.YAxis.FontSize = 24;
saveas(h6,'Change_tauF_SEIRb.png');



function dy = sens(~,y,p)
S=y(1); E=y(2); I=y(3); R=y(4); F=y(5);
dy = zeros(5,1);

% Differential equations for S, E, I, R, F
dy(1) = -p.beta/((1+p.a*F)^p.gamma)*S*I;
dy(2) =  p.beta/((1+p.a*F)^p.gamma)*S*I-E/p.tauE;
dy(3) =  E/p.tauE-I/p.tauI;
dy(4) = I/p.tauI;
dy(5) = (I-F)/p.tauF;

end

