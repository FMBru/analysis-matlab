%% gamma-2t fittness function
function out=gamma2t_minuit(par, data)
    out = par(1).*((data(1,:)-par(3))./par(2)).^2.*exp(-2.*(data(1,:)-par(3))./par(2));
    if (size(data,1)==2)    %chi-square, error = 1
        out = sum((data(2,:) - out).^2);
    elseif (size(data,1)>2) %chi-square, error = 3rd row of data
        out = sum(((data(2,:) - out)./data(3,:)).^2);
    end
end

