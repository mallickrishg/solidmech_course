% make plots of 2-d displacement fields for different displacement gradients
% Rishav Mallick, EOS, 2021

clear

x = linspace(-5,5,20);
y = linspace(-5,5,10)';

[X,Y] = meshgrid(x,y);
X0 = 0;Y0 = 0;


figure(1),clf
subplot(2,2,1)
ux = 1*(X-X0);
uy = 1*(Y-Y0);
pcolor(x,y,sqrt(ux.^2 + uy.^2)),shading interp, hold on
quiver(X,Y,ux,uy,'k','LineWidth',2)
axis tight equal, grid on
title('Extension')
colormap('gray')

subplot(2,2,2)
ux = 1*(X-X0);
uy = -1*(Y-Y0);
pcolor(x,y,sqrt(ux.^2 + uy.^2)), shading interp, hold on
quiver(X,Y,ux,uy,'k','LineWidth',2)
axis tight equal, grid on
title('Pure Shear')

subplot(2,2,3)
ux = 1*(Y-Y0);
uy = -1*(X-X0);
pcolor(x,y,sqrt(ux.^2 + uy.^2)),shading interp, hold on
quiver(X,Y,ux,uy,'k','LineWidth',2)
axis tight equal, grid on
title('Rotation')

subplot(2,2,4)
ux = 1*(Y-Y0);
uy = 1*(X-X0);
pcolor(x,y,sqrt(ux.^2 + uy.^2)),shading interp, hold on
quiver(X,Y,ux,uy,'k','LineWidth',2)
axis tight equal, grid on
title('Simple Shear')
