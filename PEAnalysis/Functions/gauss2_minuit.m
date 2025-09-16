%% gauss2 fittness function
function out=gauss2_minuit(par, data)
    %out = par(1)/par(2)/sqrt(2*pi)*exp(-1/2*((data(1,:)-par(3))/par(2)).^2) + par(4)/par(5)/sqrt(2*pi)*exp(-1/2*((data(1,:)-par(6))/par(5)).^2);
    out = par(1)*(par(4)*exp(-1/2*((data(1,:)-par(3))/par(2)).^2) + (1-par(4))*exp(-1/2*((data(1,:)-par(3))/par(5)).^2));
    if (size(data,1)==2)    %chi-square, error = 1
        out = sum((((data(2,:) - out)).^2)./out);
        %out = sum((((data(2,:) - out)).^2));
    elseif (size(data,1)>2) %chi-square, error = 3rd row of data
        out = sum(((data(2,:) - out)./data(3,:)).^2);
    end
end