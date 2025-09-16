function [xCentroid,yCentroid] = weightedAvg(Xgrid,Ygrid,weights)
%WEIGHTEDAVG: compute the weighted average on a grid given some weights

xCentroid = sum(weights(:) .* Xgrid(:)) / sum(weights(:));
yCentroid = sum(weights(:) .* Ygrid(:)) / sum(weights(:));

end

