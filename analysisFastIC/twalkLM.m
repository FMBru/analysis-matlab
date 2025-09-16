%% twalk correction function
function out = twalkLM(x, xdata)
    
    out = x(1) + x(2)./(xdata(1,:).^x(3));
    
end