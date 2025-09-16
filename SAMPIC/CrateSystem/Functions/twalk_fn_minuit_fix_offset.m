%% twalk correction function
function out = twalk_fn_minuit_fix_offset(par, fit_data)
    % Extract parameters
    % par(1) and par(2) are the variable parameters
    % fixedParam is the fixed parameter
    a = fit_data(4)
    b = par(1)
    c = par(2)
    
    % Unpack fit_data
    e_peak = fit_data(1,:);
    mean_sat = fit_data(2,:);
    e_peak_err = fit_data(3,:);
    
    % Calculate the model
    model = a + b ./ (e_peak .^ c);
    
    % Calculate the chi-square
    if (size(fit_data, 1) == 2)  % No error provided
        out = sum((mean_sat - model) .^ 2);
    elseif (size(fit_data, 1) > 2)  % Error is provided
        out = sum(((mean_sat - model) ./ e_peak_err) .^ 2);
    end
end
