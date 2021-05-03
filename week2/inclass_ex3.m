% script to compute strains from displacement field
% Rishav Mallick, EOS, 2018

clear

% load data
dstr=readtable('lateral_spread_utm.txt');
x = dstr{:,2};
y = dstr{:,3};
dux = dstr{:,4};
duy = dstr{:,6};


%% create grids and compute strain

ngridx = 50;
ngridy = 50;
nn = 10;

rthresh = 250;

xg = linspace(min(x),max(x),ngridx);
yg = linspace(min(y),max(y),ngridy);
[Xg,Yg]=meshgrid(xg,yg);

exx = zeros(numel(Xg),1);exy=exx;eyy=exx;
tic
parfor i = 1:numel(Xg)
    r2 = (Xg(i)-x).^2 + (Yg(i)-y).^2;
    [rsort,Isort]=sort(r2);
    
    if rsort(nn+1) <= (rthresh^2)       
        % create DU - 2nnx1
        Ux = dux(Isort(1)) - dux(Isort(2:nn+1));
        Uy = duy(Isort(1)) - duy(Isort(2:nn+1));
        DU = zeros(2*nn,1);
        DU(1:2:end,1) = Ux;DU(2:2:end,1) = Uy;
        % create X 2nnx4 strain components
        dx = x(Isort(1)) - x(Isort(2:nn+1));
        dy = y(Isort(1)) - y(Isort(2:nn+1));
        DX = zeros(2*nn,4);
        DX(1:2:end,1:2) = [dx dy];DX(2:2:end,3:4) = [dx dy];
        %invert matrix for displacement gradients F -
        %dux/dx,dux/dy,duy/dx,duy/dy
        F = DX\DU;
        % strain tensor exx,exy,eyy
        E = [1 0 0 0;0 0.5 0.5 0;0 0 0 1]*F;
    else
    % no strains for points with not enough neighbours    
        E = nan(3,1);
    end
    exx(i,1) = E(1);exy(i,1) = E(2);eyy(i,1) = E(3);
end
%%% plot strains
emax = sqrt(((exx-eyy)/2).^2 + exy.^2);
ep1 = (exx+eyy)/2 + emax;
ep2 = (exx+eyy)/2 - emax;
theta = 0.5*atan2d(2*exy,exx-eyy);
toc

%% plot strain components
figure(1),clf
subplot(1,4,1)
pcolor(xg./1e3,yg./1e3,reshape(ep1,ngridy,ngridx)),shading interp,axis tight equal
colormap(parula(30))
cb=colorbar;cb.Label.String = '\epsilon_1';cb.Location = 'northoutside';cb.LineWidth=1;
caxis([-1 1].*0.01)
set(gca,'FontSize',15,'LineWidth',2)

subplot(1,4,2)
pcolor(xg./1e3,yg./1e3,reshape((ep2),ngridy,ngridx)),shading interp,axis tight equal
cb=colorbar;cb.Label.String = '\epsilon_2';cb.Location = 'northoutside';cb.LineWidth=1;
caxis([-1 1].*0.01)
set(gca,'FontSize',15,'LineWidth',2)

subplot(1,4,3)
pcolor(xg./1e3,yg./1e3,reshape(emax,ngridy,ngridx)),shading interp,axis tight equal
cb=colorbar;cb.Label.String = '\epsilon_{max}';cb.Location = 'northoutside';cb.LineWidth=1;
caxis([0 1].*0.01)
set(gca,'FontSize',15,'LineWidth',2)

subplot(1,4,4)
pcolor(xg./1e3,yg./1e3,reshape(ep1+ep2,ngridy,ngridx)),shading interp,axis tight equal
cb=colorbar;cb.Label.String = '\epsilon_{kk}';cb.Location = 'northoutside';cb.LineWidth=1;
caxis([-1 1].*0.01)
set(gca,'FontSize',15,'LineWidth',2)
%% plot strain axes
box_reg = [823.5, 9890.2;...
           826.2, 9894.6];
% for plotting the strain tensor
epminval = 0.2;

figure(2),clf
pcolor(xg./1e3,yg./1e3,reshape(ep1+ep2,ngridy,ngridx)),shading interp,axis tight equal, hold on
alpha 0.5
plot_stressaxes(Xg(:),Yg(:),9,ep2,ep1,theta,1,epminval)
plot([box_reg(:,1);flipud(box_reg(:,1));box_reg(1,1)],[box_reg(1,2);box_reg(:,2);flipud(box_reg(:,2))],'k-','LineWidth',2)
cb=colorbar;cb.Label.String = '\epsilon_{kk}';cb.Location = 'northoutside';cb.LineWidth=1;
caxis([-1 1].*0.025)
xlabel('Northing (km)'), ylabel('Easting (km)')
set(gca,'FontSize',15,'LineWidth',2)

%% zoomed in plot
figure(3),clf
ind = Xg(:)./1e3>box_reg(1,1) & Xg(:)./1e3<box_reg(2,1) &...
    Yg(:)./1e3>box_reg(1,2) & Yg(:)./1e3<box_reg(2,2);

pcolor(xg./1e3,yg./1e3,reshape(ep1+ep2,ngridy,ngridx)),shading interp,axis tight equal, hold on
alpha 0.5
plot_stressaxes(Xg(ind),Yg(ind),7,ep2(ind),ep1(ind),theta(ind),0.5,epminval)
cb=colorbar;cb.Label.String = '\epsilon_{kk}';cb.Location = 'northoutside';cb.LineWidth=1;
caxis([-1 1].*0.025)
xlabel('Northing (km)'), ylabel('Easting (km)')
set(gca,'FontSize',15,'LineWidth',2)
xlim(box_reg(:,1)), ylim(box_reg(:,2))




