%% sigmoid fittness function
function out=sigmoid_fit_t(xq,yq,p)
    ycalc = ferdirac5t(xq,p);
    out = (ycalc-yq).'*(ycalc-yq); % calculate square fittnes
end