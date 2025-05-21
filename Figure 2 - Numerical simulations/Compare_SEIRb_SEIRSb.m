clearvars; close all; clc

%%%%%% Parameters and Initial Conditiosn
p.beta  = .7;
p.a = 100/3000*10^6; %adjusting for rescaling
p.gamma = 2;
p.tauE  = 5;
p.tauI  = 10;
p.tauF  = 20;
p.tauR  = 90;

S0 = 1-0.002;
E0 = 0.001;
I0 = 0.001;
R0 = 0;
F0 = 0;


IC = [S0, E0, I0, R0, F0];

%%%%%% Solvers

tmax=400000;
[t,y] = ode45(@seirb, [0:1:tmax], IC, [], p);
[t1,y1]=ode45(@seirsb,[0:1:tmax], IC, [], p);
ymax = max(max(y(:,3)),max(y1(:,3))); % Extract max value to use when setting graph limits

%%%%%%% Figures

h1=figure(1)
plot(t,y(:,3),'k',t1,y1(:,3),'--r','LineWidth',2)
xlim([0,500]);
ylim([0,ymax]);
xlabel('Time (days)','FontSize',24)
ylabel('Infectious population size (i)','FontSize',24)
legend('SEIRb','SEIRSb','FontSize',22)
legend boxoff
ax = gca;
ax.XAxis.FontSize = 24;
ax.YAxis.FontSize = 24;
print(gcf,'Compare_SEIRb_SEIRSb_500.pdf','-dpdf','-r300');
saveas(h1,'Compare_SEIRb_SEIRSb_500.png');

h2=figure(2)
plot(t,y(:,3),'k',t1,y1(:,3),'--r','LineWidth',2)
xlim([500,tmax]);
ylim([0,ymax]);
xlabel('Time (days)','FontSize',24)
xticks([500 400000])
xticklabels({'500','400,000'});
ylabel('Infectious population size (i)','FontSize',24)
legend('SEIRb','SEIRSb','FontSize',22)
legend boxoff
ax = gca;
ax.XAxis.FontSize = 24;
ax.YAxis.FontSize = 24;
print(gcf,'Compare_SEIRb_SEIRSb_equil.pdf','-dpdf','-r300');
saveas(h2,'Compare_SEIRb_SEIRSb_equil.png');


function dy = seirb(~,y,p)

s=y(1); e=y(2); i=y(3); r=y(4); f=y(5); 
dy = zeros(5,1);

% Differential equations for S, E, I, R, F
dy(1) = -p.beta/((1+p.a*f)^p.gamma)*s*i;
dy(2) =  p.beta/((1+p.a*f)^p.gamma)*s*i-e/p.tauE;
dy(3) =  e/p.tauE-i/p.tauI;
dy(4) =  i/p.tauI;
dy(5) =  (i-f)/p.tauF;

end

function dy = seirsb(~,y,p)

s=y(1); e=y(2); i=y(3); r=y(4); f=y(5); 
dy = zeros(5,1);

% Differential equations for S, E, I, R, F
dy(1) = -p.beta/((1+p.a*f)^p.gamma)*s*i+r/p.tauR;
dy(2) =  p.beta/((1+p.a*f)^p.gamma)*s*i-e/p.tauE;
dy(3) =  e/p.tauE-i/p.tauI;
dy(4) = i/p.tauI-r/p.tauR;
dy(5) = (i-f)/p.tauF;

end


