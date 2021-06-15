function [teq,Mageq,totslip,L] = compute_earthquakestatistics(ss,index,t,V,slip)
% ss - fault object
% index - fault position (in case of multiple faults, specify which fault to extract statistics)
% t - time vector
% V - velocity timeseries
% slip - cumulative slip timeseries
% Rishav Mallick, 2020, EOS

% Statistics of earthquakes
% index = ss.y2c==fltpos;

% start earthquake record from this value (yrs)
Tstart = 0;%100;%max(t)./3.15e7-900;

Vmax = max(V(:,index),[],2);
Veq = 1e-2;

teq = [];
Mageq = [];
totslip = [];
L = [];
count = 1;

FL = sum(ss.Wf);

for i = 1:length(Vmax)
    
    if Vmax(i)>Veq && Vmax(i-1)<=Veq && t(i)/3.15e7>Tstart
        teq(count,1) = t(i);
        
        slipst = slip(i,:);
        j=i;
        while Vmax(j)>Veq && j<length(Vmax)
            j = j+1;
        end
        slipend = slip(j,:);
        
        totslip(count,:) = (slipend-slipst);
        maxslip = max(totslip(count,:));
                        
        Mageq(count,1) = FL.*sum(totslip(count,:));
        
        Lindex = totslip(count,:) > 0.05*max(maxslip);
        L(count,:) = logical(Lindex);

        count = count+1;

    end
    
end


end