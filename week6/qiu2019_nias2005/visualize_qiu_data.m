clear
pc = readtable('inverse_okada_slip_dist_Nias2005_exp222_LonLatSlipTri.txt'); 
pp = readtable('niasAS_270618.txt');
ng = 300;
xg = linspace(min(pc{:,1}),max(pc{:,1}),ng);
yg = linspace(min(pc{:,2}),max(pc{:,2}),ng);
[X,Y] = meshgrid(xg,yg);
%% grid data
Fcs = scatteredInterpolant(pc{:,1},pc{:,2},pc{:,4}/100);
Fas = scatteredInterpolant(pp{:,1},pp{:,2},pp{:,3});

% bounding box
k = boundary(pc{:,1},pc{:,2});
in = inpolygon(X(:),Y(:),pc{k,1},pc{k,2});

csg = Fcs(X,Y);
csg(~in) = nan;

asg = Fas(X,Y);
asg(~in) = nan;

figure(1),clf
subplot(2,1,1)
% scatter(pc{:,1},pc{:,2},50,pp{:,3},'filled'), hold on
% plot(pc{k,1},pc{k,2},'r-','LineWidth',2)
pcolor(X,Y,csg), shading flat
axis tight equal
grid on, box on
colorbar
caxis([0 10])

subplot(2,1,2)
pcolor(X,Y,asg), shading flat
axis tight equal
grid on, box on
colorbar
caxis([0 3])
colormap(flipud(hot(10)))
%% contour it
figure(2),clf
[c_cs,hcs] = contour(X,Y,csg,[2:2:20],'LineWidth',2);
hold on
[c_as,has] = contour(X,Y,asg,[0.8:0.5:10],'LineWidth',2);

axis tight equal
grid on, box on
colorbar
caxis([0 3])
colormap(flipud(jet(15)))

% export contours
out_csfilename = 'nias_2005_cs_contours.gmt';
out_asfilename = 'nias_2005_as_contours.gmt';

% coseismic contours

fileID = fopen(out_csfilename,'w');
i = 1;
while i <=length(c_cs(1,:))
    stpt = i+1;
    npt = c_cs(2,i);
    fprintf(fileID,'%s %.2f\n','>',c_cs(1,i));
    for j = stpt:(stpt+npt-1)
        fprintf(fileID,'%.4f %.4f\n',c_cs(1,j),c_cs(2,j));
    end
    i = j+1;
end
fclose(fileID);

% afterslip contours

fileID = fopen(out_asfilename,'w');
i = 1;
while i <=length(c_as(1,:))
    stpt = i+1;
    npt = c_as(2,i);
    fprintf(fileID,'%s %.2f\n','>',c_as(1,i));
    for j = stpt:(stpt+npt-1)
        fprintf(fileID,'%.4f %.4f\n',c_as(1,j),c_as(2,j));
    end
    i = j+1;
end
fclose(fileID);
%% extract transect
addpath ~/Dropbox/scripts/utils/
w = 0.2;
x1 = 95.6;y1 = 1.7;
x2 = 97.2;y2 = 3.25;
[xp,yp] = create_swath_box(x1,y1,x2,y2,w);

figure(2),clf
[c_cs,hcs] = contour(X,Y,csg,[2:2:20],'LineWidth',2);
hold on
[c_as,has] = contour(X,Y,asg,[0.8:0.5:10],'LineWidth',2);
hold on, axis tight equal, grid on
plot([x1;x2],[y1;y2],'k-','LineWidth',2)
plot([xp;xp(1)],[yp;yp(1)],'ko-','LineWidth',1)

in = inpolygon(X(:),Y(:),[xp;xp(1)],[yp;yp(1)]);

r = 110.*sqrt((X(in)-min(X(in))).^2 + (Y(in)-min(Y(in))).^2);

figure(3),clf
plot(r,csg(in)./max(csg(in)),'k.'), hold on
plot(r,asg(in)./max(asg(in)),'b.')
axis tight

% make a transect/profile for plotting
[rsort,Isort] = sort(r);
csgsort = csg(in);csgsort = csgsort(Isort);
asgsort = asg(in);asgsort = asgsort(Isort);

nr = 100;
rvec = linspace(min(r),max(r),nr)';
dr = rvec(2)-rvec(1);
csvec = zeros(nr,1);
asvec = zeros(nr,1);

for i = 1:nr
    inr = abs(rsort-rvec(i)) <= dr;
    csvec(i) = mean(csgsort(inr));
    asvec(i) = mean(asgsort(inr));
end
plot(rvec,[csvec./max(csg(in)) asvec./max(asg(in))],'LineWidth',2)

writetable(table(rvec,[csvec./max(csg(in)) asvec./max(asg(in))]),'nias_CS_AS_transect.dat','WriteVariableNames',0,'Delimiter','\t')
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    