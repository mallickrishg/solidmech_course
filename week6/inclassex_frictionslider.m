% numerical solution for a model of linear velocity-dependent frictional
% spring-slider system
% Rishav Mallick, EOS, 2021

clear

% spring stiffness
k = 3;
% loading velocity
Vpl = 1;

% friction coeff, tau_f = A*sigma*v
Asigma = 0.5;

% relaxation time (yrs)
tR = Asigma/k;
Tend = 10*tR;%

% evl data structure
evl = [];
evl.k = k;
evl.Vpl = Vpl;
evl.Asigma = Asigma;
% set friction law
evl.frictionlaw = 1;% linear v-dependent
%% initialise vector

vi = 100*Vpl;
Y0 = vi;

% initialize the function handle 
yp=@(t,y) odefric_evo(t,y,evl);

%% Solve the system
tic
options=odeset('Refine',1,'RelTol',1e-5,'InitialStep',1e-6,'MaxStep',1e1); 
[t,Y]=ode45(yp,[0 Tend],Y0,options);
toc

% plot solution
figure(1),clf
plot(t,Y./Vpl,'x','LineWidth',2), hold on
axis tight, grid on
xlabel('t'),ylabel('v/v_{pl}')
set(gca,'Fontsize',15,'YScale','log')
plot(tR.*[1 1],get(gca,'YLim'),'k-','LineWidth',2)

