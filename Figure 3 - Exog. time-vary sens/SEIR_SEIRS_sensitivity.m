clearvars; clc; close all;

p.beta  = .7;
p.tauE  = 5;
p.tauI  = 10;
p.tauR  = 90;

s0 = (10^6-20)/10^6; 
e0 = 10/10^6; 
i0 = 10/10^6; 
r0 = 0;

IC1 = [s0, e0, i0, r0, zeros(1,4)*10e-5];
IC2 = [s0, e0, i0, r0, zeros(1,8)*10e-5];

tmax=350;
[t1,y1] = ode15s(@sens_SEIR, [0:1:tmax], IC1, 'NonNegative', p);
[t2,y2] = ode15s(@sens_SEIRS, [0:1:tmax], IC2, 'NonNegative', p);

%%%%%%% Sensitivity

h1=figure(1)
plot(t1,y1(:,7),'k',t2,y2(:,7),'r',t2,y2(:,11),'--r','LineWidth',2)
xlim([0,tmax]);
ylabel('Sensitivity','FontSize',24)
xlabel('Time (days)','FontSize',24)
legend('$i_{\beta}$','$i_{\beta}$','$i_{\tau_R}$','interpreter','latex','FontSize',28)
legend boxoff
ax = gca;
ax.XAxis.FontSize = 24;
ax.YAxis.FontSize = 24;
saveas(h1,'Sens_i_SEIR_SEIRS.png');

%%
%%%%%%%% Semi-relative (SR) sensitivity

h2=figure(2)
plot(t1,p.beta*y1(:,7),'k',t2,p.beta*y2(:,7),'r',t2,p.tauR*y2(:,11),'--r','LineWidth',2)
xlim([0,tmax]);
ylabel('Semi-relative sensitivity','FontSize',24)
xlabel('Time (days)','FontSize',24)
legend('$\beta i_{\beta}$','$\beta i_{\beta}$','$\tau_R i_{\tau_R}$','interpreter','latex','FontSize',28)
legend boxoff
ax = gca;
ax.XAxis.FontSize = 24;
ax.YAxis.FontSize = 24;
saveas(h2,'SR_Sens_i_SEIR_SEIRS.png');
%%%%%%

%% Graph all sol'ns 

h3=figure(3)
plot(t1,y1(:,3),'k',t2,y2(:,3),'r','LineWidth',2)
xlim([0,tmax]);
ylabel('Infectious population size (i)','FontSize',24)
xlabel('Time (days)','FontSize',24)
legend('SEIR','SEIRS','FontSize',22)
legend boxoff
ax = gca;
ax.XAxis.FontSize = 24;
ax.YAxis.FontSize = 24;
saveas(h3,'Plot_i_SEIR_SEIRS.png');
hold on

%%

function dy = sens_SEIR(t,y,p)
s=y(1); e=y(2); i=y(3); r=y(4);
s_beta=y(5); e_beta=y(6); i_beta=y(7); r_beta=y(8);

dy = zeros(8,1);

% Differential equations for S, E, I, R, F
dy(1) = -p.beta*s*i;
dy(2) =  p.beta*s*i-e/p.tauE;
dy(3) =  e/p.tauE-i/p.tauI;
dy(4) =  i/p.tauI;

% Partial derivatives of S, E, I, R wrt beta
dy(5) = -s*i - p.beta*i*s_beta - p.beta*s*i_beta;
dy(6) =  s*i + p.beta*i*s_beta - e_beta/p.tauE + p.beta*s*i_beta;
dy(7) =  e_beta/p.tauE - i_beta/p.tauI;
dy(8) =  i_beta/p.tauI;

end

function dy = sens_SEIRS(t,y,p)
s=y(1);         e=y(2);         i=y(3);         r=y(4);
s_beta=y(5);    e_beta=y(6);    i_beta=y(7);    r_beta=y(8);
s_tauR=y(9);    e_tauR=y(10);   i_tauR=y(11);   r_tauR=y(12);

dy = zeros(12,1);

% Differential equations for S, E, I, R
dy(1) = -p.beta*s*i+r/p.tauR;
dy(2) =  p.beta*s*i-e/p.tauE;
dy(3) =  e/p.tauE-i/p.tauI;
dy(4) =  i/p.tauI-r/p.tauR;

% Partial derivatives of S, E, I, R wrt beta
dy(5) = -s*i - p.beta*i*s_beta - p.beta*s*i_beta + r_beta/p.tauR;
dy(6) =  s*i + p.beta*i*s_beta - e_beta/p.tauE + p.beta*s*i_beta;
dy(7) =  e_beta/p.tauE - i_beta/p.tauI;
dy(8) =  i_beta/p.tauI-r_beta/p.tauR;

% Partial derivatives of S, E, I, R, F wrt beta
dy(9) = -r/p.tauR^2 - p.beta*i*s_tauR - p.beta*s*i_tauR + r_tauR/p.tauR;
dy(10) =  p.beta*i*s_tauR - e_tauR/p.tauE + p.beta*s*i_tauR;
dy(11) =  e_tauR/p.tauE - i_tauR/p.tauI;
dy(12) =  r/p.tauR^2+i_tauR/p.tauI-r_tauR/p.tauR;

end
