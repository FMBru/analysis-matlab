%% twalk correction function
function out=twalk_fn_minuit(par, data)
    
    out = par(1) + par(2)./(data(1,:).^par(3));
    if (size(data,1)==2)    %chi-square, error = 1
        out = sum((data(2,:) - out).^2);
    elseif (size(data,1)>2) %chi-square, error = 3rd row of data
        out = sum(((data(2,:) - out)./data(3,:)).^2);
    end
end