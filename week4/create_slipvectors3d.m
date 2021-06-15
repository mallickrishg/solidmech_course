% script to generate fault plane slip vectors for a known stress state
% Rishav Mallick, 2021, EOS

clear
% rng(42)
n = @(strike,dip) [cosd(strike)*sind(dip);-sind(strike)*sind(dip);cosd(dip)];
r = @(strike,dip,rake) [cosd(rake).*sind(strike)-sind(rake).*cosd(strike).*cosd(dip);...
    cosd(rake).*cosd(strike)+sind(rake).*sind(strike).*cosd(dip);...
    sind(rake).*sind(dip)]; 
s = @(strike) [sind(strike);cosd(strike);0];

% design matrix

Gfun = @(n) [n(1)-n(1)^3+n(1)*n(3)^2, n(2)-2*n(2)*n(1)^2, n(3)-2*n(3)*n(1)^2, -n(1)*n(2)^2+n(1)*n(3)^2, -2*n(1)*n(2)*n(3);...
            -n(2)*n(1)^2+n(2)*n(3)^2, n(1)-2*n(1)*n(2)^2, -2*n(1)*n(2)*n(3),   n(2)-n(2)^3+n(2)*n(3)^2,  n(3)-2*n(3)*n(2)^2;...
            -n(3)*n(1)^2-n(3)+n(3)^3, -2*n(1)*n(2)*n(3),  n(1)-2*n(1)*n(3)^2, -n(2)^2*n(3)-n(3)+n(3)^3,  n(2)-2*n(2)*n(3)^2];

        
% true stress state
s11p = -1;
s12p = -1;
s13p = -2;
s22p = 0;
s23p = 0;
s33p = -(s11p+s22p);

(s22p-s33p)/(s11p-s33p)
% sample fault planes at random
nobs = 500;

strikeobs = 180.*rand(nobs,1);
dipobs = 0 + 90.*rand(nobs,1);

%% create design matrix
G = zeros(nobs*3,5);
rhatvec = zeros(nobs*3,1);
nhatobs = zeros(nobs*3,1);
shatobs = zeros(nobs*3,1);
rakeobs = zeros(nobs,1);
rhatcheck = zeros(nobs*3,1);

for i = 1:nobs
    G(3*i-2:3*i,:) = Gfun(n(strikeobs(i),dipobs(i)));
    nhatobs(3*i-2:3*i,1) = n(strikeobs(i),dipobs(i));
    shatobs(3*i-2:3*i,1) = s(strikeobs(i));
end

% predict rake vectors
rvec = G*[s11p;s12p;s13p;s22p;s23p];
for i = 1:nobs
    rhatvec(3*i-2:3*i,1) = rvec(3*i-2:3*i,1)./sqrt(sum(rvec(3*i-2:3*i,1).^2));
    
%     rakeobs(i) = acosd(dot(rhatvec(3*i-2:3*i,1),shatobs(3*i-2:3*i,1)));
    rdots = dot(rhatvec(3*i-2:3*i,1),shatobs(3*i-2:3*i,1));
    rakeobs(i) = atan2d(rhatvec(3*i),rdots*sind(dipobs(i)));
    
    rhatcheck(3*i-2:3*i,1) = r(strikeobs(i),dipobs(i),rakeobs(i));
end

dummy = rhatvec - rhatcheck;

figure(11),clf
plot3(rhatvec(1:3:end),rhatvec(2:3:end),rhatvec(3:3:end),'kp'), hold on
plot(cosd(linspace(0,180)),sind(linspace(0,180)),'k-','LineWidth',2)
plot(cosd(linspace(0,180)),-sind(linspace(0,180)),'k-','LineWidth',2)
axis tight equal
grid on, box on
view(0,90)

T = [zeros(nobs,2),strikeobs,dipobs,rakeobs];
save('INPUTdata.mat','T')