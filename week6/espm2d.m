function [ux,uz] = espm2d(m,x)
% calculate the velocity field due to interseismic locking using an elastic
% subducting plate model (Kanda and Simons, 2010)
% INPUTS 
% m - vector containing [Vpl,dip(º),D(km),T(km)]
% x - station distance from the trench (km)
% OUTPUTS
% ux - horizontal velocity (trench perpendicular)
% uz - vertical velocity 
% Rishav Mallick, EOS 2018

zeta = @(x,xd,D) (x-xd)./D;

uzdd = @(m,x) m(1)/pi.*(sind(m(2)).*atan2(zeta(x,m(3)/tand(m(2)),m(3)),1) + (cosd(m(2))+zeta(x,m(3)/tand(m(2)),m(3)).*sind(m(2)))./(1+zeta(x,m(3)/tand(m(2)),m(3)).^2));
uzdslab = @(m,x) -m(1)/pi.*(sind(m(2)).*atan2(zeta(x,-m(4)*tand(m(2)/2),m(4)),1) + (cosd(m(2))+zeta(x,-m(4)*tand(m(2)/2),m(4)).*sind(m(2)))./(1+zeta(x,-m(4)*tand(m(2)/2),m(4)).^2));
uzfslab = @(m,x) -m(1)/pi.*(-1./(1+zeta(x,-m(4)*tand(m(2)/2),m(4)).^2));

uxdd = @(m,x) m(1)/pi.*(cosd(m(2)).*atan2(zeta(x,m(3)/tand(m(2)),m(3)),1) + (sind(m(2))-zeta(x,m(3)/tand(m(2)),m(3)).*cosd(m(2)))./(1+zeta(x,m(3)/tand(m(2)),m(3)).^2));
uxdslab = @(m,x) -m(1)/pi.*(cosd(m(2)).*atan2(zeta(x,-m(4)*tand(m(2)/2),m(4)),1) + (sind(m(2))-zeta(x,-m(4)*tand(m(2)/2),m(4)).*cosd(m(2)))./(1+zeta(x,-m(4)*tand(m(2)/2),m(4)).^2));
uxfslab = @(m,x) -m(1)/pi.*(-atan2(zeta(x,-m(4)*tand(m(2)/2),m(4)),1) + (zeta(x,-m(4)*tand(m(2)/2),m(4)))./(1+zeta(x,-m(4)*tand(m(2)/2),m(4)).^2));

uzESPM = @(m,x) uzdd(m,x) + uzdslab(m,x) + uzfslab(m,x);
uxESPM = @(m,x) uxdd(m,x) + uxdslab(m,x) + uxfslab(m,x);

uz = uzESPM(m,x);
ux = uxESPM(m,x);
end