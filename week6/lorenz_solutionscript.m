% lorenz attractor ODE solution
clear

evl = [];
evl.sigma = 10;
evl.b = 8/3;
evl.r = 28;

%% initial conditions

Y0 = [0,1,3];

% initialize the function handle 
yp=@(t,y) odelorenz(t,y,evl);

Tend = 50;
% tspan = [0:0.001:50];
%% numerical solution
tic
options=odeset('Refine',1,'RelTol',1e-12,'AbsTol',1e-12*ones(1,3)); 
[t,Y]=ode45(yp,[0 Tend],Y0,options);
% [t,Y]=ode45(yp,tspan,Y0,options);
toc

xvec = Y(:,1);
yvec = Y(:,2);
zvec = Y(:,3);

%% plot solution
figure(1),clf
plot(t,xvec,'LineWidth',1), hold on
plot(t,yvec+40,'LineWidth',1)
plot(t,zvec-40,'LineWidth',1)
legend('x','y','z')
axis tight
% xlim([0 Tend])
xlabel('Time')
set(gca,'FontSize',15)

%% plot butterfly's wings
figure(10),clf
plot3(xvec,yvec,zvec,'w-','LineWidth',1)
axis tight equal, box on, grid on
xlabel('x'),ylabel('y'),zlabel('z')
set(gca,'FontSize',20,'LineWidth',2,'Color','k','GridColor','w')   
view(47,19)
return
%% plot as movie
figure(2),clf
ni = 2000;
for i = 1+ni:40:length(xvec)
    plot3(xvec(i-ni:i),yvec(i-ni:i),zvec(i-ni:i),'k-','LineWidth',1)
    axis tight, box on, grid on
    
    xlim([min(xvec) max(xvec)])
    ylim([min(yvec) max(yvec)])
    zlim([min(zvec) max(zvec)])
    
    pause(0.01)
end