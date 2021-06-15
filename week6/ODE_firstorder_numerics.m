% script to demonstrate different numerical integration methods for ODEs
% Rishav Mallick, 2021, EOS

clear

% input parameters for ODE
Vpl = 1;
k = 3;
tauinf = k*Vpl;
Asigma = 0.5;

tr = Asigma/k;

dydt = @(t,v) (tauinf - k*v)/(Asigma);

% initial condition
vi = 10;

% Analytical solution
ytrue = @(t,vi) vi*exp(-t./tr) + Vpl*(1-exp(-t./tr));
%% solve using numerical methods

nsteps = 5;
t = linspace(0,5*tr,nsteps)';
dt = [0;diff(t)];

%% Euler's method

yeuler = zeros(nsteps,1);

yeuler(1) = vi;

for j=2:nsteps
    
    h = dt(j);
    
    yeuler(j,1) = yeuler(j-1,1) + h*dydt(t(j-1,1),yeuler(j-1,1));
    
end

%% Heun's Method

yheun = zeros(nsteps,1);
yheun(1) = vi;

for j=2:nsteps
    
    h = dt(j);
    
    k1 = h*dydt( t(j-1), yheun(j-1) );
    
    k2 = h*dydt( t(j-1)+h, yheun(j-1)+k1);
    
    yheun(j,1) = yheun(j-1) + 0.5*(k1+k2);
end

%% RK 4th order method

yrk = zeros(nsteps,1);
yrk(1) = vi;


for j=2:nsteps
    
    h = dt(j);
    
    k1 = h*dydt( t(j-1), yrk(j-1) );
    
    k2 = h*dydt( t(j-1)+h/2, yrk(j-1)+0.5*k1 );
    
    k3 = h*dydt( t(j-1)+h/2, yrk(j-1)+0.5*k2 );
    
    k4 = h*dydt( t(j-1)+h, yrk(j-1)+k3 );
    
    yrk(j,1) = yrk(j-1)+(1/6)*(k1+2*k2+2*k3+k4);
end

%% plot solutions

figure(1),clf
plot(t,yeuler,'-','LineWidth',2,'Color',rgb('steelblue')), hold on
plot(t,yheun,'-','LineWidth',2,'Color',rgb('peru'))
plot(t,yrk,'-','LineWidth',2,'Color',rgb('forestgreen'))

plot(t,ytrue(t,vi),'r--','LineWidth',1)
plot([1 1].*tr,get(gca,'YLim'),'k-','LineWidth',2)

axis tight, grid on
xlabel('t'),ylabel('v(t)')
set(gca,'FontSize',15,'YScale','lin','LineWidth',1)













