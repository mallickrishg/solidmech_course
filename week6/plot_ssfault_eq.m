% plot surface displacements for a strike slip earthquake

%  demo on how to use unicycle
% Rishav Mallick, EOS, 2021

clear
addpath ~/Dropbox/scripts/utils/
addpath ../unicycle/matlab/
import unicycle.*

% read fault geometry file
G = 30e3;
nu = 0.25;
patchfname = 'ssflt_3D.seg';
earthModel = unicycle.greens.okada92(G,nu);
rcv = unicycle.geometry.receiver(patchfname,earthModel);

dipdist = cumsum(rcv.W,1);

%% calculate dispalcement kernels
nx = 102;
ny = 100;

ox = linspace(-250e3,250e3,nx)';
oy = linspace(-250e3,250e3,ny)';
[X,Y] = meshgrid(ox,oy);

% displacement kernels
[Gs,Gd] = rcv.displacementKernels([X(:),Y(:),0.*X(:)],3);

slip = zeros(rcv.N,1);
meanc = mean(rcv.xc);
dist_meanc = (rcv.xc(:,1)-meanc(1)).^2 + (rcv.xc(:,2)-meanc(2)).^2 + (rcv.xc(:,3)-meanc(3)).^2;
ind = dist_meanc==min(dist_meanc);
slip(ind) = 1;

% displacements
u = Gs*slip;
ue = u(1:3:end);
un = u(2:3:end);
uz = u(3:3:end);

% downsample vectors
nq = 13;

%% plot surface displacements

figure(1),clf

subplot(121)
imagesc(ox./1e3,oy./1e3,reshape(uz,ny,nx))
hold on
contour(ox./1e3,oy./1e3,reshape(uz,ny,nx),[-1:0.1:1].*1e-2,'k')
quiver(X(1:nq:end)'./1e3,Y(1:nq:end)'./1e3,ue(1:nq:end),un(1:nq:end),'k-','LineWidth',2)
axis tight equal, grid on, box on
colormap jet(200)
cb=colorbar;cb.Label.String = 'u_z';
caxis([-1 1].*1e-2)
xlabel('x_2 (km)'),ylabel('x_1 (km)')
set(gca,'YDir','normal','FontSize',20)

subplot(122)
imagesc(ox./1e3,oy./1e3,reshape(uz,ny,nx))
hold on
contour(ox./1e3,oy./1e3,reshape(uz,ny,nx),[-1:0.1:1].*5e-4,'k')
quiver(X(1:nq:end)'./1e3,Y(1:nq:end)'./1e3,ue(1:nq:end),un(1:nq:end),'k-','LineWidth',2)
axis tight equal, grid on, box on
colormap jet(200)
cb=colorbar;cb.Label.String = 'u_z';
caxis([-1 1].*5e-4)
xlabel('x_2 (km)'),ylabel('x_1 (km)')
set(gca,'YDir','normal','FontSize',20)

print('ssfault_surfdisplacements','-djpeg','-r300')