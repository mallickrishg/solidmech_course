clear

load eulercalc_matlab/countries_map.mat
load eulercalc_matlab/GPS_ref.mat
%% extract data in box
box_lon = [74 95];
box_lat = [25 38];

% project data from spherical coordinates to local x-y coordinate system
[x, y] = latlon_to_xy_polyconic(lat, lon, mean(box_lat), mean(box_lon));
[xf, yf] = latlon_to_xy_polyconic(Flat, Flon, mean(box_lat), mean(box_lon));
[xcoast,ycoast] = latlon_to_xy_polyconic(coastlines(:,2), coastlines(:,1), mean(box_lat), mean(box_lon));

in_box = lon<=box_lon(2) & lon>=box_lon(1) & lat>=box_lat(1) & lat<=box_lat(2);
flt_box = (Flon<=box_lon(2) & Flon>=box_lon(1) & Flat>=box_lat(1) & Flat<=box_lat(2)) | isnan(Flon);
coast_box = coastlines(:,1)<=box_lon(2) & coastlines(:,1)>=box_lon(1) & coastlines(:,2)>=box_lat(1) & coastlines(:,2)<=box_lat(2);

figure(1),clf
plot(coastlines(:,1),coastlines(:,2),'k-','LineWidth',1), hold on
plot(Flon,Flat,'r-')
quiver(lon(in_box),lat(in_box),uE(in_box),uN(in_box),'LineWidth',2)
axis tight equal, grid on
xlim(box_lon)
ylim(box_lat)
xlabel('Lon º'), ylabel('Lat º')
set(gca,'FontSize',15,'LineWidth',2)

figure(2),clf
plot(xcoast(coast_box),ycoast(coast_box),'k.','LineWidth',1), hold on
plot(xf(flt_box),yf(flt_box),'r-')
quiver(x(in_box),y(in_box),uE(in_box),uN(in_box),'LineWidth',2)
axis tight equal, grid on
xlabel('x (km)'),ylabel('y (km)')
set(gca,'FontSize',15,'LineWidth',2)

%% rotate coordinates
rot = -20;
[x_r,y_r] = rotate_xy(x(in_box),y(in_box),rot);
[ux_r,uy_r] = rotate_xy(uE(in_box),uN(in_box),rot);
[xf_r,yf_r] = rotate_xy(xf(flt_box),yf(flt_box),rot);

figure(3),clf
plot(xf_r,yf_r,'r-','LineWidth',1), hold on
quiver(x_r,y_r,ux_r,uy_r,'LineWidth',2)
axis tight equal, box on, grid on
xlabel('x (km)'),ylabel('y (km)')
set(gca,'FontSize',15,'LineWidth',2)

