function Yp = odefric_evo(~,Y,evl)

if evl.frictionlaw == 1
    v = Y;
    mv = (evl.m*v^(1-1/evl.m));
    %mv = evl.m*(v^0.5);
    Yp = mv*evl.k*(evl.Vpl-v)./(evl.Asigma);
end

end