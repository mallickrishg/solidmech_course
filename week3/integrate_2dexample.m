% homework assignment - integration of 1-d and 2-d strain fields
% compute numerical integration and compare with analytical solution
% Rishav Mallick, EOS, 2021

clear

x = linspace(-5,5,200);
y = linspace(0.5,10,200);
[X,Y] = meshgrid(x,y);

% strain functions
ezx = X./(X.^2 + Y.^2);
ezy = Y./(X.^2 + Y.^2);

% extract left boundary
leftbound = X==min(x);
%extract bottom boundary
botbound = Y==min(y);

uznumleft = zeros(size(X));
uznumleft(leftbound) = cumtrapz(Y(leftbound),2.*ezy(leftbound));

% calculate numerical solution
uznum = zeros(size(X));
% use uznumleft as a boundary condition
for i = 1:length(y)
    uznum(i,:) = cumtrapz(x,2.*ezx(i,:)) + uznumleft(i,1);
end

% analytical solution
uz = log(X.^2 + Y.^2);

% plot strains
figure(2),clf
subplot(121)
imagesc(x,y,ezx), hold on
contourf(x,y,ezx,20,'k-')
axis tight equal
title('\epsilon_{xz}');
cb=colorbar;
xlabel('x'),ylabel('y')
set(gca,'YDir','normal','FontSize',20)

subplot(122)
imagesc(x,y,ezy), hold on
contourf(x,y,ezy,20,'k-')
title('\epsilon_{yz}');
axis tight equal
cb=colorbar;
xlabel('x'),ylabel('y')
set(gca,'YDir','normal','FontSize',20)

% plot displacement and comparison
figure(3),clf
subplot(121)
imagesc(x,y,uz), hold on
contourf(x,y,uz,20,'k-')
axis tight equal
cb=colorbar;
colormap(flipud(hot(200)))
title('u_z (analytical)')
xlabel('x'), ylabel('y')
set(gca,'YDir','normal','FontSize',20)

subplot(122)
imagesc(x,y,uznum), hold on
contourf(x,y,uznum,20,'k-')
axis tight equal
cb=colorbar;
title('u_z (numerical)')
set(gca,'YDir','normal','FontSize',20)


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

