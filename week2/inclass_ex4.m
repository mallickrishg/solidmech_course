% generate dispalcement fields from potential function and use it to
% visualize strains
% Rishav Mallick, EOS, 2021

clear

% potential functions
phi = @(m,a,x2,x3) 0.5*m.*log(x2.^2 + x3.^2) + 2.*sin(a*x3.*x2);
u2_a = @(m,a,x2,x3) m.*x2./(x2.^2 + x3.^2) + 2*a*x3.*cos(a.*x3.*x2);
u3_a = @(m,a,x2,x3) m.*x3./(x2.^2 + x3.^2) + 2*a*x2.*cos(a.*x3.*x2);

m = 3;
a = 0;

x2 = linspace(-3,3,90)';
x3 = linspace(-3,3,100);

[X2,X3] = meshgrid(x2,x3);

potfield = phi(m,a,X2,X3);
[u2,u3] = gradient(potfield,x2,x3);
%% plot potential field and displacements
np = 3;
nsamp = 500;
rng(42)

figure(1),clf
imagesc(x2,x3,potfield), hold on
% contour(x2,x3,potfield,10,'k-','LineWidth',1)
startx2 = -3 + 6.*rand(nsamp,1);
startx3 = -3 + 6.*rand(nsamp,1);

hl=streamline(X2,X3,u2,u3,startx2,startx3);
set(hl,'Linewidth',1,'Color','w')
quiver(X2(1:np:end),X3(1:np:end),u2(1:np:end),u3(1:np:end),'r','LineWidth',1)
axis tight equal
caxis([-1 1].*5)
cb=colorbar;cb.Label.String = '\phi(x_2,x_3)';
set(gca,'FontSize',15,'LineWidth',1)


figure(2),clf
subplot(2,3,1)
imagesc(x2,x3,u2), hold on
contour(x2,x3,u2,[-5:5],'k','LineWidth',1)
axis tight equal
cb=colorbar;cb.Label.String = 'u_2';
set(gca,'FontSize',15,'LineWidth',1)
caxis([-1 1].*5)

subplot(2,3,4)
imagesc(x2,x3,u3), hold on
contour(x2,x3,u3,[-5:5],'k','LineWidth',1)
axis tight equal
cb=colorbar;cb.Label.String = 'u_3';
set(gca,'FontSize',15,'LineWidth',1)
caxis([-1 1].*5)

subplot(2,3,2)
imagesc(x2,x3,u2_a(m,a,X2,X3)), hold on
contour(x2,x3,u2_a(m,a,X2,X3),[-5:5],'k','LineWidth',1)
axis tight equal
cb=colorbar;cb.Label.String = 'u_2 (analytical)';
set(gca,'FontSize',15,'LineWidth',1)
caxis([-1 1].*5)

subplot(2,3,5)
imagesc(x2,x3,u3_a(m,a,X2,X3)), hold on
contour(x2,x3,u3_a(m,a,X2,X3),[-5:5],'k','LineWidth',1)
axis tight equal
cb=colorbar;cb.Label.String = 'u_3 (analytical)';
set(gca,'FontSize',15,'LineWidth',1)
caxis([-1 1].*5)

subplot(2,3,3)
imagesc(x2,x3,u2_a(m,a,X2,X3)-u2)
axis tight equal
cb=colorbar;cb.Label.String = 'error (u_2)';
set(gca,'FontSize',15,'LineWidth',1)
caxis([-1 1].*1)

subplot(2,3,6)
imagesc(x2,x3,u3_a(m,a,X2,X3)-u3)
axis tight equal
cb=colorbar;cb.Label.String = 'error (u_3)';
set(gca,'FontSize',15,'LineWidth',1)
caxis([-1 1].*1)

%% compute strains
% [e22,u23] = gradient(u2_a(m,a,X2,X3),x2,x3);
% [u32,e33] = gradient(u3_a(m,a,X2,X3),x2,x3);
[e22,u23] = gradient(u2,x2,x3);
[u32,e33] = gradient(u3,x2,x3);

e23 = 0.5*(u23 + u32);
curlu = u23-u32;
divu = e22+e33;

emax = sqrt(((e22-e33)/2).^2 + e23.^2);
ep1 = (e22+e33)/2 + emax;
ep2 = (e22+e33)/2 - emax;
theta = 0.5*atan2d(2*e23,e22-e33);
epminval = 0.05;

figure(3),clf
imagesc(x2,x3,ep1+ep2)
alpha 0.5
plot_stressaxes(X2.*1e3,X3.*1e3,6,ep2,ep1,theta,1,epminval)
axis tight equal
cb=colorbar;cb.Label.String = '\epsilon_{kk}';
set(gca,'FontSize',15,'LineWidth',1)
caxis([-1 1].*10)

figure(4),clf
subplot(121)
imagesc(x2,x3,divu), hold on
contour(x2,x3,divu,[-10:2:10],'k-','LineWidth',1)
quiver(X2(1:np:end),X3(1:np:end),u2(1:np:end),u3(1:np:end),'r','LineWidth',1)
axis tight equal
cb=colorbar;cb.Label.String = '\nabla.u';
set(gca,'FontSize',15,'LineWidth',1)
caxis([-1 1].*10)

subplot(122)
imagesc(x2,x3,curlu), hold on
quiver(X2(1:np:end),X3(1:np:end),u2(1:np:end),u3(1:np:end),'r','LineWidth',1)
axis tight equal
cb=colorbar;cb.Label.String = '\nabla x u';
set(gca,'FontSize',15,'LineWidth',1)
caxis([-1 1].*1)








