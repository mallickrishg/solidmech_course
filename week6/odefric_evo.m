function Yp = odefric_evo(~,Y,evl)

if evl.frictionlaw == 1
    v = Y;
    Yp = evl.k*(evl.Vpl-v)./(evl.Asigma);
end

end