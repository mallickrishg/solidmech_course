function [vE,vN] = pole_velocity(lat,lon,latp,lonp,degmyr)
%function to calculate the predicted velocities at lat,lon (sites) given the euler
%pole - latp,lonp,degmyr
%Geodetic (ellipsoid) coordinates in degrees are assumed and converted to
%cartesian coordinates using a best fit spherical geometry for the ellipsoid. 
% Rishav Mallick, 2017, EOS

%define spheroid constraints from WGS1984
fi=1/298.257223563;
e2=2*fi-fi^2;
% convert geodetic latitudes to bestfit sphere
lat=atan((1-e2)*tan(pi/180*(lat)))*180/pi;
latp=atan((1-e2)*tan(pi/180*(latp)))*180/pi;
% do rotation
omega=euler_vector(latp,lonp,degmyr);
Rx=rot_matrix(lat,lon);
V=Rx*omega;
vE=V(1:2:end);
vN=V(2:2:end);

end