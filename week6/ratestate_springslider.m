% numerical solution for a model of rate-state-dependent frictional
% spring-slider system
% Rishav Mallick, EOS, 2021
clc
clear

% evl data structure
evl = [];
evl.K = 1;% MPa/m
evl.sigma = 10;% MPa

evl.Vpl = 1e-8;
evl.Vo = 1e-6;

% friction properties
evl.a = 1e-2;
evl.b = evl.a + 0.5e-2;
evl.l = 1e-2;
evl.mu0 = 0.6;

% critical spring stiffness
hstar = (evl.b-evl.a)*evl.sigma/evl.l;

disp(['k_c = ' num2str(hstar)])
disp(['k = ' num2str(evl.K)])
%% initialise vector

s_i = 0;
psi_i = 0;
zeta_i = -4;

Y0 = [s_i,psi_i,zeta_i];

% initialize the function handle 
yp=@(t,y) ode_rsf(t,y,evl);

%% numerical solution

Tend = 10*3.15e7;

tic
options=odeset('Refine',1,'RelTol',1e-6,'MaxStep',1e8); 
[t,Y]=ode45(yp,[0 Tend],Y0,options);
toc

%% extract outputs
slip = Y(:,1);
f = evl.mu0 + evl.a*Y(:,3) + evl.b*Y(:,2);
V = evl.Vo*exp(Y(:,3));

figure(1),clf
subplot(3,1,1)
plot(t./3.15e7,slip,'-','LineWidth',2)
axis tight, grid on
set(gca,'FontSize',15,'LineWidth',1)
xlabel('Time (yr)'), ylabel('Slip (m)')

subplot(312)
semilogy(t./3.15e7,V,'-','LineWidth',2)
axis tight, grid on
set(gca,'FontSize',15,'LineWidth',1)
xlabel('Time (yr)'), ylabel('V (m/s)')

subplot(313)
plot(t./3.15e7,f,'-','LineWidth',2)
axis tight, grid on
set(gca,'FontSize',15,'LineWidth',1)
xlabel('Time (yr)'), ylabel('f')


figure(2),clf
subplot(211)
plot(V,f,'LineWidth',2)
axis tight, grid on
ylabel('f'), xlabel('v')
ylim([0.5 0.7])
set(gca,'FontSize',15,'LineWidth',1,'XScale','log')

subplot(212)
plot(exp(Y(1,2))*evl.l/evl.Vo,V(1),'kp','MarkerSize',10), hold on
plot(exp(Y(:,2))*evl.l/evl.Vo,V,'LineWidth',2)
axis tight, grid on
xlabel('\theta'), ylabel('v')
set(gca,'FontSize',15,'LineWidth',1,'XScale','log','YScale','log')


figure(3),clf
semilogy(slip,V,'k-','LineWidth',2)
axis tight, grid on
ylabel('V')
yyaxis right
plot(slip,f,'LineWidth',2),axis tight
ylabel('f')
xlim([1.2 2.2])
xlabel('Slip')
set(gca,'FontSize',20)










