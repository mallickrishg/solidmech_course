function plot_stressaxes(X,Z,np,Sigmin,Sigmax,theta,scf,Rminval)
% plot stress axes given the following quantities
% Sigmin,Sigmax,theta (principal stress field orientation and magnitude)
% X,Z - gridded values
% scf - scaling factor for plotting axes
% Rminval - select a threshold value for Sigmin or Sigmax (for plotting)
% Rishav Mallick, EOS, 2020

if scf==0
    Rmin = Sigmin(1:np:end,1:np:end);
    Rmin = abs(Rmin)./Rmin;
    Rmax = Sigmax(1:np:end,1:np:end);
    Rmax = abs(Rmax)./Rmax;
    qscale = 0.2;
else
    Rmin = Sigmin(1:np:end,1:np:end).*scf;
    Rmin = Rmin./max(abs(Rmin(:)));
    Rmax = Sigmax(1:np:end,1:np:end).*scf;
    Rmax = Rmax./max(abs(Rmax(:)));
    qscale = 0.4*scf;
end

Rmin = Rmin(:);Rmax = Rmax(:);
thetaplot = theta(1:np:end,1:np:end);thetaplot = thetaplot(:);

xplot = X(1:np:end,1:np:end)./1e3; xplot = xplot(:);
zplot = Z(1:np:end,1:np:end)./1e3; zplot = zplot(:);

% filter out small values
% Rminval = 0.1;
plotid = abs(Rmax)<Rminval & abs(Rmin)<Rminval;
xplot(plotid) = [];
zplot(plotid) = [];
Rmin(plotid) = [];
Rmax(plotid) = [];
thetaplot(plotid) = [];

% remove nans
ind = isnan(thetaplot);

% condition checking for orientation of strain tensor
IF_extcomp = Rmax(:)>0 & Rmin(:)<0; %Sigmax extensional & Sigmin is compressional
IF_extext = Rmax(:)>0 & Rmin(:)>0;%Sigmax extensional & Sigmin is extensional
IF_compcomp = Rmax(:)<0 & Rmin(:)<0; %Sigmax compressional & Sigmin is compressional

% plot compressional stresses (Sigmax and Sigmin)
l = 2;
w = 0.3;

%% extension_compression
xp = xplot(IF_extcomp & ~ind);
zp = zplot(IF_extcomp & ~ind);
Rpc = Rmax(IF_extcomp & ~ind).*cosd(thetaplot(IF_extcomp & ~ind));
Rps = Rmax(IF_extcomp & ~ind).*sind(thetaplot(IF_extcomp & ~ind));
arrow([xp,zp],[xp+qscale.*Rpc,zp+qscale.*Rps],'Length',l,'Width',w)
arrow([xp,zp],[xp-qscale.*Rpc,zp-qscale.*Rps],'Length',l,'Width',w)
Rpc = Rmin(IF_extcomp & ~ind).*cosd(thetaplot(IF_extcomp & ~ind)+90);
Rps = Rmin(IF_extcomp & ~ind).*sind(thetaplot(IF_extcomp & ~ind)+90);
arrow([xp+qscale.*Rpc,zp+qscale.*Rps],[xp,zp],'Length',l,'Width',w)
arrow([xp-qscale.*Rpc,zp-qscale.*Rps],[xp,zp],'Length',l,'Width',w)

%% extension_extension
xp = xplot(IF_extext & ~ind);
zp = zplot(IF_extext & ~ind);
Rpc = Rmax(IF_extext & ~ind).*cosd(thetaplot(IF_extext & ~ind));
Rps = Rmax(IF_extext & ~ind).*sind(thetaplot(IF_extext & ~ind));
arrow([xp,zp],[xp+qscale.*Rpc,zp+qscale.*Rps],'Length',l,'Width',w)
arrow([xp,zp],[xp-qscale.*Rpc,zp-qscale.*Rps],'Length',l,'Width',w)
Rpc = Rmin(IF_extext & ~ind).*cosd(thetaplot(IF_extext & ~ind)+90);
Rps = Rmin(IF_extext & ~ind).*sind(thetaplot(IF_extext & ~ind)+90);
arrow([xp,zp],[xp+qscale.*Rpc,zp+qscale.*Rps],'Length',l,'Width',w)
arrow([xp,zp],[xp-qscale.*Rpc,zp-qscale.*Rps],'Length',l,'Width',w)

%% compression_compression

xp = xplot(IF_compcomp & ~ind);
zp = zplot(IF_compcomp & ~ind);
Rpc = Rmax(IF_compcomp & ~ind).*cosd(thetaplot(IF_compcomp & ~ind));
Rps = Rmax(IF_compcomp & ~ind).*sind(thetaplot(IF_compcomp & ~ind));
arrow([xp+qscale.*Rpc,zp+qscale.*Rps],[xp,zp],'Length',l,'Width',w)
arrow([xp-qscale.*Rpc,zp-qscale.*Rps],[xp,zp],'Length',l,'Width',w)
Rpc = Rmin(IF_compcomp & ~ind).*cosd(thetaplot(IF_compcomp & ~ind)+90);
Rps = Rmin(IF_compcomp & ~ind).*sind(thetaplot(IF_compcomp & ~ind)+90);
arrow([xp+qscale.*Rpc,zp+qscale.*Rps],[xp,zp],'Length',l,'Width',w)
arrow([xp-qscale.*Rpc,zp-qscale.*Rps],[xp,zp],'Length',l,'Width',w)


end