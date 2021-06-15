%  demo on how to use unicycle
% Rishav Mallick, EOS, 2021

clear
addpath ~/Dropbox/scripts/utils/
addpath ../unicycle/matlab/
import unicycle.*

% read fault geometry file
G = 30e3;
nu = 0.25;
patchfname = 'flt_3D.seg';
earthModel = unicycle.greens.okada92(G,nu);
rcv = unicycle.geometry.receiver(patchfname,earthModel);

dipdist = cumsum(rcv.W,1);
%% compute stresses and displacements using Okada1992

k = 5;
n2 = 240;
n3 = 70;

x2 = linspace(-20e3,100e3,n2);
x3 = linspace(-35e3,-0.1,n3);
[X2,X3] = meshgrid(x2,x3);

% compute distance from top left corner
% xd=rcv.xc(:,1)-rcv.x(k,1);
% yd=rcv.xc(:,2)-rcv.x(k,2);
% zd = rcv.xc(:,3);
xd=X2(:)-rcv.x(k,1);
yd=X2(:).*0 -rcv.x(k,2);
% yd=rcv.L(1)/2 + zeros(numel(X2),1);
zd=X3(:);    

M = numel(X2);

[U,~,~,S] = unicycle.greens.computeOkada92(1,xd(:),yd(:),zd(:),G,nu, ...
        -rcv.x(k,3),rcv.dip(k)/180*pi,rcv.L(k),rcv.W(k),'d',0,rcv.strike(k)/180*pi);
S=reshape(S,M,3,3);
U=reshape(U,M,3);

% compute displacement components
u1 = U(:,1);
u2 = U(:,2);
u3 = U(:,3);

% compute stress components
S11 = S(:,1,1);
S12 = S(:,1,2);
S13 = S(:,1,3);
S22 = S(:,2,2);
S23 = S(:,2,3);
S33 = S(:,3,3);

figure(1),clf
imagesc(x2./1e3,x3./1e3,reshape(S11,n3,n2)), hold on
plot(rcv.xc(:,1)./1e3,rcv.xc(:,3)./1e3,'k-')
% quiver(X2(:)./1e3,X3(:)./1e3,u1,u3,'k-','LineWidth',2)
axis tight equal
colormap bluewhitered
cb=colorbar;
xlabel('x_2 (km)'),ylabel('x_3 (km)')
set(gca,'YDir','normal','FontSize',20)

%% compute traction kernels
[Kss,Kds,Ksd,Kdd,~,~]=rcv.tractionKernels(rcv);

%% surface displacements
ox = linspace(30e3,100e3,1000)';

[Gs,Gd] = rcv.displacementKernels([ox,0.*ox,0.*ox],3);

Ge = Gd(1:3:end,:);
Gn = Gs(2:3:end,:);
Gz = Gd(3:3:end,:);



































%% export slip distribution
slip = zeros(rcv.N,1);
taudrop = zeros(rcv.N,1);
in = dipdist>=10e3 & dipdist<=80e3;
taudrop(in) = -2;

slip(in) = Kdd(in,in)\taudrop(in);

uz = Gz*slip;
ue = Ge*slip;

figure(11),clf
plot(ox./1e3,uz,'-','LineWidth',2), hold on
plot(ox./1e3,ue,'-','LineWidth',2)
% plot(ox./1e3,Gn*slip,'k-','LineWidth',2)
axis tight, box on, grid on

T = table(ox,ue,uz);
% writetable(T,'eq_planestrain.dat')

