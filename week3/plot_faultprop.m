function plot_faultprop(flt,prop)

if length(prop(1,:)) == 1
    
    toplot = prop;
    
    for i = 1:flt.N
        patch([flt.xTL(i,1),flt.xTR(i,1),flt.xBR(i,1),flt.xBL(i,1),flt.xTL(i,1)]./1e3,[flt.xTL(i,2),flt.xTR(i,2),flt.xBR(i,2),flt.xBL(i,2),flt.xTL(i,2)]./1e3,...
            [flt.xTL(i,3),flt.xTR(i,3),flt.xBR(i,3),flt.xBL(i,3),flt.xTL(i,3)]./1e3,toplot(i)), hold on
    end
    
else
    
    toplot = sqrt(prop(:,1).^2 + prop(:,2).^2);
    ss1 = prop(:,2);
    ss2 = prop(:,2);
    for i = 1:flt.N
        patch([flt.xTL(i,1),flt.xTR(i,1),flt.xBR(i,1),flt.xBL(i,1),flt.xTL(i,1)]./1e3,[flt.xTL(i,2),flt.xTR(i,2),flt.xBR(i,2),flt.xBL(i,2),flt.xTL(i,2)]./1e3,...
            [flt.xTL(i,3),flt.xTR(i,3),flt.xBR(i,3),flt.xBL(i,3),flt.xTL(i,3)]./1e3,toplot(i)), hold on
    end
    
    scf = 1;
    
    slipvec = flt.sv.*repmat(ss1,1,3) - flt.dv.*repmat(ss2,1,3);
    % plot slipvectors
    quiver3(flt.xc(:,1)./1e3,flt.xc(:,2)./1e3,flt.xc(:,3)./1e3,slipvec(:,1).*scf,slipvec(:,2).*scf,slipvec(:,3).*scf,'b-','LineWidth',2)
end

colormap(hot(100))
end