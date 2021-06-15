% script to demonstrate Taylor series expansion
% Rishav Mallick, EOS, 2021

clear

ng = 100;
x = linspace(-1*pi,1*pi,ng)';

% specify function to be approximated
y = @(x) sin(x);
% about which point
x0 = pi/2;

% Taylor series
taylor1 = y(x0) + (x-x0).*cos(x0);
taylor2 = y(x0) + (x-x0).*cos(x0) - sin(x0)/factorial(2).*(x-x0).^2;
taylor3 = y(x0) + (x-x0).*cos(x0) - sin(x0)/factorial(2).*(x-x0).^2 - cos(x0)/factorial(3).*(x-x0).^3;
taylor4 = y(x0) + (x-x0).*cos(x0) - sin(x0)/factorial(2).*(x-x0).^2 - cos(x0)/factorial(3).*(x-x0).^3 + sin(x0)/factorial(4).*(x-x0).^4;
taylor5 = y(x0) + (x-x0).*cos(x0) - sin(x0)/factorial(2).*(x-x0).^2 - cos(x0)/factorial(3).*(x-x0).^3 + sin(x0)/factorial(4).*(x-x0).^4 + cos(x0)/factorial(5).*(x-x0).^5;

taylorexp = y(x0) + (x-x0).*cos(x0) - sin(x0)/factorial(2).*(x-x0).^2 ...
     - cos(x0)/factorial(3).*(x-x0).^3 + sin(x0)/factorial(4).*(x-x0).^4 ...
    + cos(x0)/factorial(5).*(x-x0).^5 - sin(x0)/factorial(6).*(x-x0).^6 ...
    - cos(x0)/factorial(7).*(x-x0).^7 + sin(x0)/factorial(8).*(x-x0).^8 ...
    + cos(x0)/factorial(9).*((x-x0).^9) - sin(x0)./factorial(10).*(x-x0).^10;

%% plot figure comparing different approximations

figure(1),clf

plot(x,y(x),'k-','LineWidth',3), hold on
plot(x0,y(x0),'rp','MarkerFaceColor','r','MarkerSize',20)
% plot(x,taylor1,'-','LineWidth',2)
% plot(x,taylor2,'-','LineWidth',2)
% plot(x,taylor3,'-','LineWidth',2)
plot(x,taylor4,'-','LineWidth',2)
% plot(x,taylor5,'-','LineWidth',2)

% plot(x,taylorexp,'--','LineWidth',3)

axis tight, grid on
xlabel('x'),ylabel('y(x)')
set(gca,'FontSize',25)