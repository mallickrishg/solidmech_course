% Script for in-class exercise 
% Rishav Mallick, EOS, 2021
clear
%% calculating numerical derivatives

% in 1-D
x = linspace(-7,7,100)';
f = sin(x) + log(x.^2);

figure(1),clf
plot(x,f,'k-','LineWidth',2)
axis tight, grid on, box on
xlabel('x'),ylabel('f(x)')
set(gca,'Fontsize',15)

% dydx = [nan;diff(f)./diff(x)];
dydx = gradient(f,x);
% analytical derivative
dydx_a = cos(x) + 2./x;

yyaxis right
plot(x,dydx,'-','LineWidth',2), hold on
plot(x,dydx_a,'r--','LineWidth',2)
axis tight
xlabel('x'),ylabel('df/dx')

figure(2),clf
subplot(122)
plot(x(x>0),dydx(x>0),'-','LineWidth',2), hold on
plot(x(x>0),dydx_a(x>0),'r--','LineWidth',2)
plot(x(x>0),0.*x(x>0),'k-')
axis tight, grid on
xlabel('x'),ylabel('df/dx')
ylim([-2 max(abs(dydx))])
set(gca,'Fontsize',15)

subplot(121)
plot(x(x<=0),dydx(x<=0),'-','LineWidth',2), hold on
plot(x(x<=0),dydx_a(x<=0),'r--','LineWidth',2)
plot(x(x<=0),0.*x(x<=0),'k-')
axis tight, grid on
xlabel('x'),ylabel('df/dx')
ylim([-max(abs(dydx)) 2])
set(gca,'Fontsize',15)

%% 2-D using finite differences

x = linspace(-7,7,50)';
y = linspace(-4,4,50)';
[X,Y] = meshgrid(x,y);

f = sin(X.*Y) + log(X.^2);

figure(10),clf
imagesc(x,y,f)
axis tight equal
cb = colorbar;cb.Label.String = 'f(x,y)';
xlabel('x'),ylabel('y')
set(gca,'FontSize',15)

% numerical gradient
[dfdx,dfdy]=gradient(f,x,y);
%analytical gradient
dfdx_a = Y.*cos(X.*Y) + 2./X;
dfdy_a = X.*cos(X.*Y);

figure(11),clf
% numerical
subplot(231)
imagesc(x,y,dfdx)
axis tight equal
cb = colorbar;cb.Label.String = 'f_x(x,y)';
xlabel('x'),ylabel('y')
set(gca,'FontSize',15)
title('Numerical gradient')

subplot(234)
imagesc(x,y,dfdy)
axis tight equal
cb = colorbar;cb.Label.String = 'f_y(x,y)';
xlabel('x'),ylabel('y')
set(gca,'FontSize',15)

% analytical
subplot(232)
imagesc(x,y,dfdx_a)
axis tight equal
cb = colorbar;cb.Label.String = 'f_x(x,y)';
xlabel('x'),ylabel('y')
set(gca,'FontSize',15)
title('Analytical gradient')

subplot(235)
imagesc(x,y,dfdy_a)
axis tight equal
cb = colorbar;cb.Label.String = 'f_y(x,y)';
xlabel('x'),ylabel('y')
set(gca,'FontSize',15)

% error (analytical - numerical)
subplot(233)
imagesc(x,y,dfdx_a-dfdx)
axis tight equal
cb = colorbar;cb.Label.String = 'error_x';
xlabel('x'),ylabel('y')
set(gca,'FontSize',15)
title('Analytical - Numerical')

subplot(236)
imagesc(x,y,dfdy_a-dfdy)
axis tight equal
cb = colorbar;cb.Label.String = 'error_y';
xlabel('x'),ylabel('y')
set(gca,'FontSize',15)




