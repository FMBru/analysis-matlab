function [meanAmp,meanSAT, stdSAT, meanriseTime, stdriseTime, twalk_par, tres_par, riseTime_par, riseTimeStd_par] = characterisePad(table,padNumber, folder, shouldSave, opt)
%CHARACTERISEPAD 

if nargin < 5 || isempty(opt)
    opt = 'pad';
end

vars = {'MM_amp', 'SAT', 'e_peak_MM', 'riseTime'};
padNumberString = num2str(padNumber);

if strcmp(opt, 'full')
    for i=1:length(vars)
        vars{i} = [vars{i} padNumberString];
    end
elseif strcmp(opt, 'central')
    padNumberString = [padNumberString '-CENTRAL'];
end

padNumberString = ['-pad' padNumberString];

meanAmp = ampDistribution(table.MCP_amp,table.(vars{1}), padNumberString, folder, shouldSave);
[meanSAT, stdSAT] = distribution(table.(vars{2}),100,'SAT','ns', folder, padNumberString, shouldSave); %SAT distribution over the pad
[meanriseTime, stdriseTime] = distribution(table.(vars{4}),100,'riseTime','ns', folder, padNumberString, shouldSave);     %riseTime distribution over the pad   
[~, twalk_par, tres_par] = correlationPlots(table.(vars{3}),table.(vars{2}), 'e-peak-charge', 'pC', 'SAT', 'ns', folder, padNumberString, shouldSave, true); %SAT vs e peak over the pad
[~, riseTime_par, riseTimeStd_par] = correlationPlots(table.(vars{3}),table.(vars{4}), 'e-peak-charge', 'pC', 'riseTime', 'ns', folder, padNumberString, shouldSave, true); %riseTime vs e peak over the pad

end

