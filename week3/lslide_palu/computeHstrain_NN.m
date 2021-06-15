function [exx,eyy,exy] = computeHstrain_NN(x,y,Ue,Un,ngridx,ngridy,nn)
% compute horizontal strain tensor from horizontal displacements in
% cartesian coordinates using nearest neighbors (method is ideal for
% closely sampled noisy data)
% INPUTS -
% x,y,Ue,Un - vectors of horizontal position(x,y) and horizontal
% displacements/veclocities (Ue,Un)
% ngridx,ngridy - scalar values for number of grid points in x,y direction
% at which we will evaluate strain
% nn - number of nearby stations to use to calculate strain (greater value
% implies more smoothing)
% OUTPUTS - exx,eyy,exy - strain tensor
% Rishav Mallick, EOS, December 2018

xg = linspace(min(x),max(x),ngridx);
yg = linspace(min(y),max(y),ngridy);
[Xg,Yg]=meshgrid(xg,yg);

exx = zeros(numel(Xg),1);exy=exx;eyy=exx;

for i = 1:numel(Xg)
    r2 = (Xg(i)-x).^2 + (Yg(i)-y).^2;
    [~,Isort]=sort(r2);
    % create DU - 2nnx1
    Ux = Ue(Isort(1)) - Ue(Isort(2:nn+1));
    Uy = Un(Isort(1)) - Un(Isort(2:nn+1));
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
    exx(i,1) = E(1);exy(i,1) = E(2);eyy(i,1) = E(3);
end

end