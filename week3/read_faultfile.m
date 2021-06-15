function flt = read_faultfile(filename)

datatab = readtable(filename);

% store data in a structure
flt.xc = datatab{:,1:3};
flt.L = datatab.L;
flt.W = datatab.W;
flt.strike = datatab.strike;
flt.dip = datatab.dip;

flt.N = length(flt.L);

%% unit vectors
flt.nv = [cosd(flt.strike).*sind(flt.dip),-sind(flt.strike).*sind(flt.dip),cosd(flt.dip)];
flt.sv = [sind(flt.strike),cosd(flt.strike),zeros(flt.N,1)];
flt.dv = [cosd(flt.strike).*cosd(flt.dip),-sind(flt.strike).*cosd(flt.dip),-sind(flt.dip)];

%% plotting
flt.xTL = flt.xc - 0.5.*repmat(flt.L,1,3).*flt.sv - 0.5.*repmat(flt.W,1,3).*flt.dv;
flt.xTR = flt.xc + 0.5.*repmat(flt.L,1,3).*flt.sv - 0.5.*repmat(flt.W,1,3).*flt.dv;
flt.xBL = flt.xc - 0.5.*repmat(flt.L,1,3).*flt.sv + 0.5.*repmat(flt.W,1,3).*flt.dv;
flt.xBR = flt.xc + 0.5.*repmat(flt.L,1,3).*flt.sv + 0.5.*repmat(flt.W,1,3).*flt.dv;

end
