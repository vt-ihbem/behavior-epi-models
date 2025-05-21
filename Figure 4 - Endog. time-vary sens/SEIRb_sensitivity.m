clearvars; close all; clc

p.beta  = .7;
p.a = 100/3000*10^6; %adjusting for rescaling;
p.gamma = 2;
p.tauE  = 5;
p.tauI  = 10;
p.tauF  = 20;

s0 = (10^6-20)/10^6; 
e0 = 10/10^6; 
i0 = 10/10^6; 
r0 = 0;
f0 = 0;

IC = [s0, e0, i0, r0, f0, zeros(1,15)*10e-5];

tmax=350;
[t,y] = ode15s(@sens, [0:1:tmax], IC, [], p);

%%%%%%% Sensitivity

h1=figure(1)
plot(t,y(:,8),'k',t,y(:,13),'--k',t,y(:,18),'-.k','LineWidth',2)
xlim([0,tmax]);
ylabel('Sensitivity','FontSize',24)
xlabel('Time (days)','FontSize',24)
legend('$i_{\beta}$','$i_a$',...
        '$i_{\gamma}$','interpreter','latex','FontSize',22)
legend boxoff
ax = gca;
ax.XAxis.FontSize = 24;
ax.YAxis.FontSize = 24;
saveas(h1,'Sens_i_SEIRb.png');

%%
%%%%%%%% Semi-relative (SR) sensitivity

h2=figure(2)
plot(t,p.beta*y(:,8),'k',t,p.a*y(:,13),'--k',t,p.gamma*y(:,18),'-.k','LineWidth',2)
xlim([0,tmax]);
ylabel('Semi-relative sensitivity','FontSize',24)
xlabel('Time (days)','FontSize',24)
legend('$\beta i_{\beta}$','$a i_a$',...
        '$\gamma i_{\gamma}$','interpreter','latex','FontSize',22) %'NumColumns',2
legend boxoff
ax = gca;
ax.XAxis.FontSize = 24;
ax.YAxis.FontSize = 24;
saveas(h2,'SR_Sens_i_SEIRb.png');
%%%%%%

%% Graph all sol'ns 

h3=figure(3)
plot(t,y(:,3),'k','LineWidth',2)
xlim([0,tmax]);
ylabel('Infectious population size (i)','FontSize',24)
xlabel('Time (days)','FontSize',24)
ax = gca;
ax.XAxis.FontSize = 24;
ax.YAxis.FontSize = 24;
saveas(h3,'Plot_i_SEIRb.png');

%%

function dy = sens(t,y,p)
s=y(1); e=y(2); i=y(3); r=y(4); f=y(5); 
s_beta=y(6); e_beta=y(7); i_beta=y(8); r_beta=y(9); f_beta=y(10);
s_a=y(11); e_a=y(12); i_a=y(13); r_a=y(14); f_a=y(15);
s_gamma=y(16); e_gamma=y(17); i_gamma=y(18); r_gamma=y(19); f_gamma=y(20);
dy = zeros(20,1);

% Differential equations for S, E, I, R, F
dy(1) = -p.beta/((1+p.a*f)^p.gamma)*s*i;
dy(2) =  p.beta/((1+p.a*f)^p.gamma)*s*i-e/p.tauE;
dy(3) =  e/p.tauE-i/p.tauI;
dy(4) = i/p.tauI;
dy(5) = (i-f)/p.tauF;

% Partial derivatives of S, E, I, R, F wrt beta
dy(6) = -1/((1+p.a*f)^p.gamma)*s*i-p.beta/((1+p.a*f)^p.gamma)*i*s_beta-p.beta/((1+p.a*f)^p.gamma)*s*i_beta...
            +(p.a*p.beta*p.gamma)/((1+p.a*f)^(p.gamma+1))*s*i*f_beta;
dy(7) =  1/((1+p.a*f)^p.gamma)*s*i+p.beta/((1+p.a*f)^p.gamma)*i*s_beta-e_beta/p.tauE...
            +p.beta/((1+p.a*f)^p.gamma)*s*i_beta-(p.a*p.beta*p.gamma)/((1+p.a*f)^(p.gamma+1))*s*i*f_beta;
dy(8) = e_beta/p.tauE-i_beta/p.tauI;
dy(9) = i_beta/p.tauI;
dy(10)= (i_beta-f_beta)/p.tauF;

% Partial derivatives of S, E, I, R, F wrt a
dy(11) =  (p.beta*p.gamma*f)/((1+p.a*f)^(p.gamma+1))*s*i-p.beta/((1+p.a*f)^(p.gamma))*i*s_a...
            -p.beta/((1+p.a*f)^(p.gamma))*s*i_a+(p.a*p.beta*p.gamma)/((1+p.a*f)^(p.gamma+1))*s*i*f_a;
dy(12) = -(p.beta*p.gamma*f)/((1+p.a*f)^(p.gamma+1))*s*i+p.beta/((1+p.a*f)^p.gamma)*i*s_a-(1/p.tauE)*e_a...
              +p.beta/((1+p.a*f)^p.gamma)*s*i_a-(p.a*p.beta*p.gamma)/((1+p.a*f)^(p.gamma+1))*s*i*f_a;
dy(13) = e_a/p.tauE-i_a/p.tauI;
dy(14) = i_a/p.tauI;
dy(15)= (i_a-f_a)/p.tauF;

% Partial derivatives of S, E, I, R, F wrt gamma
dy(16) =  p.beta/((1+p.a*f)^(p.gamma))*s*i*log(1+p.a*f)-p.beta/((1+p.a*f)^(p.gamma))*i*s_gamma...
            -p.beta/((1+p.a*f)^(p.gamma))*s*i_gamma+(p.a*p.beta*p.gamma)/((1+p.a*f)^(p.gamma+1))*s*i*f_gamma;
dy(17) = -p.beta/((1+p.a*f)^(p.gamma))*s*i*log(1+p.a*f)+p.beta/((1+p.a*f)^(p.gamma))*i*s_gamma...
            -e_gamma/p.tauE+p.beta/((1+p.a*f)^(p.gamma))*s*i_gamma-(p.a*p.beta*p.gamma)/((1+p.a*f)^(p.gamma+1))*s*i*f_gamma;
dy(18) = e_gamma/p.tauE-i_gamma/p.tauI;
dy(19) = i_gamma/p.tauI;
dy(20)= (i_gamma-f_gamma)/p.tauF;
end

