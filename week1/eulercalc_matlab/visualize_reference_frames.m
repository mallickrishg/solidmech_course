% Script to visualize data in your favourite reference frame
% Provide path to the dat file and provide the euler pole you want to
% transform it to
% The data is in ITRF08 reference frame
% Rishav Mallick, Nov 2017, Myanmar workshop
close all
clear

%load international boundaries
load ../countries_map.mat
%load GPS data in ITRF08 (from Kreemer et al., 2014)
load ../GPS_ref.mat

inp_transform=input('Do you want to transform to another reference frame? (yes=1,no=0) \n');
%Ind/ITRF08 is lat=51N,lon=12E,degmyr=0.5

if inp_transform==1
    latp=input('Enter the euler pole latitude (º) ');
    lonp=input('Enter the euler pole longitude (º) ');
    degmyr=input('Enter the rotation rate (º/Myr) ');
    [uEn,uNn]=pole_velocity(lat,lon,latp,lonp,degmyr);
    figure(1),clf
    plot(coastlines(:,1),coastlines(:,2),'k'),hold on
    quiver(lon,lat,uE-uEn,uN-uNn,1)
    title('Global map of station velocities','FontSize',15)
    xlabel('Longitude (\circ)')
    ylabel('Latitude (\circ)')
    axis tight, axis equal
else
    disp('Visualizing default velocities')
    figure(1),clf
    plot(coastlines(:,1),coastlines(:,2),'k'),hold on
    quiver(lon,lat,uE,uN,1)
    title('Global map of station velocities','FontSize',15)
    xlabel('Longitude (\circ)')
    ylabel('Latitude (\circ)')
    axis tight, axis equal
end

%% zoomed in plots
XL = [70 110];
YL = [-10 50];
scf = 20;

figure(2),clf
subplot(121)
plot(coastlines(:,1),coastlines(:,2),'k'),hold on
quiver(lon,lat,uE./scf,uN./scf,0,'LineWidth',1)
% only for velocity scale bar
quiver(75,0,50./scf,0./scf,0,'LineWidth',1)
axis tight equal
set(gca,'Fontsize',15)
xlabel('lat'),ylabel('lon')
xlim(XL)
ylim(YL)
title('ITRF vel field')

subplot(122)
plot(coastlines(:,1),coastlines(:,2),'k'),hold on
quiver(lon,lat,uEn./scf,uNn./scf,0,'LineWidth',1)
% only for velocity scale bar
quiver(75,0,50./scf,0./scf,0,'LineWidth',1)
axis tight equal
set(gca,'Fontsize',15)
xlim(XL)
ylim(YL)
title('predicted motion for Eurasia Euler pole')
