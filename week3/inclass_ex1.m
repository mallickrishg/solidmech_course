% homework assignment - integration of 1-d and 2-d strain fields
% compute numerical integration and compare with analytical solution
% Rishav Mallick, EOS, 2021

clear

% 1-d function integration
x = linspace(1e-14,5,1000)';

% function def
f = log(x);

% analytical solution (without integration constant)
F = x.*(log(x)-1);

% numerical integration
Fnum = cumtrapz(x,f);

figure(1),clf
subplot(211)
plot(x,f,'LineWidth',2)
axis tight
xlabel('x')
ylabel('f(x)')
set(gca,'Fontsize',20)

subplot(212)
plot(x,F,'b-','LineWidth',2), hold on
plot(x,Fnum,'r-','LineWidth',2)
axis tight
xlabel('x')
ylabel('F = \int f dx')
set(gca,'Fontsize',20)

%% 2-d function
clear

x = linspace(-5,5,100);
y = linspace(0.5,10,90);
[X,Y] = meshgrid(x,y);

% strain functions
ezx = X./(X.^2 + Y.^2);
ezy = Y./(X.^2 + Y.^2);

% extract left boundary
leftbound = X==min(x);
%extract bottom boundary
botbound = Y==min(y);

uznumleft = zeros(size(X));
uznumleft(leftbound) = cumtrapz(Y(leftbound),ezy(leftbound));

% calculate numerical solution
uznum = zeros(size(X));
% use uznumleft as a boundary condition
for i = 1:length(y)
    uznum(i,:) = cumtrapz(x,ezx(i,:)) + uznumleft(i,1);
end

% analytical solution
uz = 1.*log(X.^2 + Y.^2);% - 1.6;

% plot strains
figure(2),clf
subplot(121)
imagesc(x,y,ezx), hold on
contour(x,y,ezx,20,'k-')
axis tight equal
title('\epsilon_{xz}');
cb=colorbar;
xlabel('x'),ylabel('y')
set(gca,'YDir','normal','FontSize',20)

subplot(122)
imagesc(x,y,ezy), hold on
contour(x,y,ezy,20,'k-')
title('\epsilon_{yz}');
axis tight equal
cb=colorbar;
xlabel('x'),ylabel('y')
set(gca,'YDir','normal','FontSize',20)

% plot displacement and comparison
figure(3),clf
subplot(121)
imagesc(x,y,uz), hold on
contour(x,y,uz,20,'k-')
axis tight equal
cb=colorbar;
colormap(flipud(hot(200)))
title('u_z (analytical)')
xlabel('x'), ylabel('y')
set(gca,'YDir','normal','FontSize',20)

subplot(122)
imagesc(x,y,uznum), hold on
contour(x,y,uznum,20,'k-')
axis tight equal
cb=colorbar;
title('u_z (numerical)')
set(gca,'YDir','normal','FontSize',20)

% plot difference between numerical and analytical solutions
figure(4),clf
imagesc(x,y,uz-uznum), hold on
contour(x,y,uz-uznum,5,'k-')
axis tight equal
cb=colorbar;colormap gray
set(gca,'YDir','normal','FontSize',15)

% plot profiles of dispalcement along boundaries
figure(5),clf
subplot(121)
plot(x,uz(botbound),'LineWidth',2), hold on
plot(x,uznum(botbound),'LineWidth',2)
axis tight, grid on
legend('analytical','numerical')
xlabel('x'),ylabel('u_z = \int\epsilon_{zx}dx')
set(gca,'YDir','normal','FontSize',20)

subplot(122)
plot(uz(leftbound),y,'LineWidth',2), hold on
plot(uznum(leftbound),y,'LineWidth',2)
axis tight, grid on
ylabel('y'),xlabel('u_z = \int\epsilon_{zy}dy')
set(gca,'FontSize',20)

%% plane-strain problem

clear

x = linspace(-5,5,200);
y = linspace(0.5,5,110);
[X,Y] = meshgrid(x,y);

k = 4;

exx = ((k+1)*(X.^2.*Y) + (k-3).*Y.^3)./(X.^2 + Y.^2).^2;
eyy = ((k-3)*(X.^2.*Y) + (k+1).*Y.^3)./(X.^2 + Y.^2).^2;
exy = 4*(X.*Y.^2)./(X.^2 + Y.^2).^2;

% analytical solutions
ux = (k-1)*atan2(X,Y) - 2*(X.*Y)./(X.^2 + Y.^2);
uy = (X.^2-Y.^2)./(X.^2 + Y.^2) + (k+1)/2.*log(X.^2 + Y.^2);

% extract left boundary
leftbound = X==min(x);
%extract bottom boundary
botbound = Y==min(y);

uxnumleft = zeros(size(X));
uynumbot = zeros(size(X));

% uxnumleft(leftbound) = cumtrapz(Y(leftbound),exy(leftbound));
uxnumleft(leftbound) = ux(leftbound);
% uynumbot(botbound) = cumtrapz(X(botbound),exy(botbound));
uynumbot(botbound) = uy(botbound);

% calculate numerical solution
uynum = zeros(size(X));
uxnum = zeros(size(X));
% use uznumleft as a boundary condition
for i = 1:length(y)
    uxnum(i,:) = cumtrapz(x,exx(i,:)) + uxnumleft(i,1);
end
for i = 1:length(x)
    uynum(:,i) = cumtrapz(y,eyy(:,i)) + uynumbot(1,i);
end


% numerical gradients
[uxx,uxy] = gradient(ux,x,y);
[uyx,uyy] = gradient(uy,x,y);

% figure(100),clf
% subplot(121)
% plot(y,uxnumleft(leftbound),'LineWidth',2), hold on
% plot(y,cumtrapz(Y(leftbound),exy(leftbound)),'r-','Linewidth',2)
% plot(y,cumtrapz(Y(leftbound),uxy(leftbound)),'b--','Linewidth',1)
% xlabel('y'), ylabel('u_x')
% axis tight, grid on
% set(gca,'FontSize',15)
% 
% subplot(122)
% plot(x,uynumbot(botbound),'LineWidth',2), hold on
% plot(x,cumtrapz(X(botbound),2*exy(botbound)),'r-','Linewidth',2)
% xlabel('x'), ylabel('u_y')
% axis tight, grid on
% set(gca,'FontSize',15)
%
% figure(22),clf
% subplot(311)
% imagesc(x,y,uxx), hold on
% contour(x,y,uxx,20,'k-')
% axis tight equal
% cb=colorbar;
% colormap(parula(200))
% title('\epsilon_{xx}')
% set(gca,'YDir','normal','FontSize',20)
% 
% subplot(312)
% imagesc(x,y,uyy), hold on
% contour(x,y,uyy,20,'k-')
% axis tight equal
% cb=colorbar;
% title('\epsilon_{yy}')
% set(gca,'YDir','normal','FontSize',20)
% 
% subplot(313)
% imagesc(x,y,(uxy+uyx)./2), hold on
% contour(x,y,(uxy+uyx)./2,20,'k-')
% axis tight equal
% cb=colorbar;
% title('\epsilon_{xy}')
% set(gca,'YDir','normal','FontSize',20)
% %
% strains
figure(10),clf
subplot(311)
imagesc(x,y,exx), hold on
contour(x,y,exx,20,'k-')
axis tight equal
cb=colorbar;
colormap(parula(200))
title('\epsilon_{xx}')
set(gca,'YDir','normal','FontSize',20)

subplot(312)
imagesc(x,y,eyy), hold on
contour(x,y,eyy,20,'k-')
axis tight equal
cb=colorbar;
title('\epsilon_{yy}')
set(gca,'YDir','normal','FontSize',20)

subplot(313)
imagesc(x,y,exy), hold on
contour(x,y,exy,20,'k-')
axis tight equal
cb=colorbar;
title('\epsilon_{xy}')
set(gca,'YDir','normal','FontSize',20)

% % numerical displacements
% figure(11),clf
% subplot(122)
% imagesc(x,y,uynum), hold on
% contour(x,y,uynum,20,'k-')
% axis tight equal
% cb=colorbar;
% colormap(hot(200))
% title('u_y (numerical)')
% set(gca,'YDir','normal','FontSize',20)
% 
% subplot(121)
% imagesc(x,y,uxnum), hold on
% contour(x,y,uxnum,20,'k-')
% axis tight equal
% cb=colorbar;
% title('u_x (numerical)')
% set(gca,'YDir','normal','FontSize',20)

% verify integration

figure(12),clf
subplot(221)
imagesc(x,y,uxnum), hold on
contour(x,y,uxnum,20,'k-')
colorbar
axis tight equal
title('u_x (numerical)')
set(gca,'YDir','normal','FontSize',20)

subplot(222)
imagesc(x,y,ux), hold on
contour(x,y,ux,20,'k-')
colorbar
axis tight equal
set(gca,'YDir','normal','FontSize',20)

subplot(223)
imagesc(x,y,uynum), hold on
contour(x,y,uynum,20,'k-')
colorbar
axis tight equal
set(gca,'YDir','normal','FontSize',20)
title('u_y (numerical)')

subplot(224)
imagesc(x,y,uy), hold on
contour(x,y,uy,20,'k-')
colorbar
axis tight equal
set(gca,'YDir','normal','FontSize',20)
colormap hot
