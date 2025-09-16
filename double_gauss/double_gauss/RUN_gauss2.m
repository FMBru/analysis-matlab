clear all;
close all;

% Define the list of filenames to analyze
filenames = {
% Photocathode studies
"Run314-MM_res_20MOhm_10_mm-120_um-PC_CsI_5_nm_REVERSED-Flushing-A275V_C455V-MCP_cut_150_mV_and_geo_cut-SATlist.txt",
"Run317-MM_res_20MOhm_10_mm-120_um-PC_CsI_5_nm_REVERSED-Flushing-A275V_C445V-MCP_cut_150_mV_and_geo_cut-SATlist.txt",
"Run318-MM_res_20MOhm_10_mm-120_um-PC_CsI_5_nm_REVERSED-Flushing-A275V_C435V-MCP_cut_150_mV_and_geo_cut-SATlist.txt",
"Run319-MM_res_20MOhm_10_mm-120_um-PC_CsI_5_nm_REVERSED-Flushing-A265V_C465V-MCP_cut_150_mV_and_geo_cut-SATlist.txt"
};

% Create a list of full file paths
file_list = cellfun(@(filename) fullfile('Input', filename), filenames, 'UniformOutput', false);
      
 
file_list_ROOT = {
% % Photocathodes studies
% fullfile('ResultsROOT', 'Run271_ROOTfit.txt')
% fullfile('ResultsROOT', 'Run273_ROOTfit.txt')
fullfile('ResultsROOT', 'Run314_ROOTfit.txt')
fullfile('ResultsROOT', 'Run317_ROOTfit.txt')
fullfile('ResultsROOT', 'Run318_ROOTfit.txt')
fullfile('ResultsROOT', 'Run319_ROOTfit.txt')
% fullfile('ResultsROOT', 'Run355_ROOTfit.txt')
% fullfile('ResultsROOT', 'Run356_ROOTfit.txt')
% fullfile('ResultsROOT', 'Run357_ROOTfit.txt')
% fullfile('ResultsROOT', 'Run359_ROOTfit.txt')
% fullfile('ResultsROOT', 'Run360_ROOTfit.txt')
% fullfile('ResultsROOT', 'Run382_ROOTfit.txt')
% fullfile('ResultsROOT', 'Run383_ROOTfit.txt')
% fullfile('ResultsROOT', 'Run390_ROOTfit.txt')
% fullfile('ResultsROOT', 'Run391_ROOTfit.txt')
% fullfile('ResultsROOT', 'Run392_ROOTfit.txt')
% fullfile('ResultsROOT', 'Run396_ROOTfit.txt')
% fullfile('ResultsROOT', 'Run398_ROOTfit.txt')
% fullfile('ResultsROOT', 'Run399_ROOTfit.txt')
% fullfile('ResultsROOT', 'Run215_ROOTfit.txt')
% fullfile('ResultsROOT', 'Run216_ROOTfit.txt')
% fullfile('ResultsROOT', 'Run217_ROOTfit.txt')
% fullfile('ResultsROOT', 'Run248_ROOTfit.txt')
% fullfile('ResultsROOT', 'Run264_ROOTfit.txt')
% fullfile('ResultsROOT', 'Run231_ROOTfit.txt')
% fullfile('ResultsROOT', 'Run232_ROOTfit.txt')
% fullfile('ResultsROOT', 'Run233_ROOTfit.txt')
% fullfile('ResultsROOT', 'Run234_ROOTfit.txt')
};

% Define the number of bins for histogram
nbins = 100;
min = -0.2;
max = 0.2;
    
    % Loop over each file in the list
for file_index = 1:numel(file_list)
    % Load the data from the current file
    data = load(file_list{file_index});
    data = data(data > min & data < max);
    
    % Read the content of the text file as a string
    fid = fopen(file_list_ROOT{file_index}, 'r');
    file_content = fscanf(fid, '%c', inf);
    fclose(fid);

    % Extract key-value pairs using regular expressions
    key_value_pairs = regexp(file_content, 'fit\.(\w+)\s*=\s*([-0-9\.]+);', 'tokens');

    
    % Initialize the struct
    fit = struct();

    % Assign values to the struct fields
    for i = 1:numel(key_value_pairs)
        key = key_value_pairs{i}{1};
        value = str2double(key_value_pairs{i}{2});
        fit.(key) = value;
    end

%     if file_index == 2
%         fit.mean = -2.7455e-05;
%     end
% 
%     if file_index == 9
%         fit.mean = -8.27408e-05;
%     end
% % 
%     if file_index == 23
%         fit.mean = -1.35153e-05;
%     end

    % Plot dual Gauss histogram
    figure;
    binss = linspace(min,max,nbins+1);
    %h = histogram(data, nbins, 'FaceColor', 'b');
    h = histogram(data, binss, 'FaceColor', [0 0.4470 0.7410]);
    hold on;
    fit_data = [h.BinEdges(1:end-1) + h.BinWidth/2; h.Values];

    ax = gca;
    mycolors = [1 0 0; 0 0.25 1; 0 0.85 0];
    ax.ColorOrder = mycolors;
    
    
    if isfield(fit, 'mean')
        plot(fit_data(1,:), gauss1_minuit([fit.scale*fit.mixing fit.sigma_core fit.mean], fit_data(1,:)), 'LineWidth', 1.5);
        plot(fit_data(1,:), gauss1_minuit([fit.scale*(1-fit.mixing) fit.sigma_tail fit.mean], fit_data(1,:)), 'LineWidth', 1.5);
        plot(fit_data(1,:), gauss2_minuit([fit.scale fit.sigma_core fit.mean fit.mixing fit.sigma_tail], fit_data(1,:)), 'LineWidth', 2);
    else    
        plot(fit_data(1,:), gauss1_minuit([fit.scale*fit.mixing fit.sigma_core fit.g1_mean], fit_data(1,:)), 'LineWidth', 1.5);
        plot(fit_data(1,:), gauss1_minuit([fit.scale*(1-fit.mixing) fit.sigma_tail fit.g1_mean], fit_data(1,:)), 'LineWidth', 1.5);
        plot(fit_data(1,:), gauss2_minuit([fit.scale fit.sigma_core fit.g1_mean fit.mixing fit.sigma_tail], fit_data(1,:)), 'LineWidth', 2);
    end
 

    message = sprintf('   Entries = %d\n',fit.Entries);
    message = [message sprintf('   \\chi^2 / NDF = %2.2f / %2d\n',fit.ChiSquare, fit.NDF)];
    if isfield(fit, 'mean')
        message = [message sprintf('   \\mu = %2.2f ps \\pm %2.2f ps\n',1000*fit.mean,1000*fit.mean_error)];
    else
        message = [message sprintf('   \\mu = %2.2f ps \\pm %2.2f ps\n',1000*fit.g1_mean,1000*fit.g1_mean_error)];
    end
    message = [message sprintf('   \\sigma_{1} = %2.2f ps \\pm %2.2f ps\n',1000*fit.sigma_core,1000*fit.sigma_core_error)];
    message = [message sprintf('   \\sigma_{2} = %2.2f ps \\pm %2.2f ps\n',1000*fit.sigma_tail,1000*fit.sigma_tail_error)];
    message = [message sprintf('   \\sigma_{comb} = %2.2f ps \\pm %2.2f ps \n',1000*fit.sigma_all,1000*fit.sigma_all_error)];
    message = [message sprintf('   RMS_{tot} = %2.2f ps \\pm %2.2f ps', 1000*fit.StdDev, 1000*fit.StdDev/sqrt(2*fit.Entries))];


    % Extract the mean value from the struct
%     mean_value = fit.g1_mean;
% 
%     % Plot the Gaussians using the extracted mean value
%     plot(fit_data(1,:), gauss1_minuit([fit.scale*fit.mixing fit.sigma_core mean_value], fit_data(1,:)), 'LineWidth', 1.5);
%     plot(fit_data(1,:), gauss1_minuit([fit.scale*(1-fit.mixing) fit.sigma_tail mean_value], fit_data(1,:)), 'LineWidth', 1.5);
%     plot(fit_data(1,:), gauss2_minuit([fit.scale fit.sigma_core mean_value fit.mixing fit.sigma_tail], fit_data(1,:)), 'LineWidth', 2);
% 
%     message = sprintf('   Entries = %d\n',fit.Entries);
%     message = [message sprintf('   \\chi^2 / NDF = %2.2f / %2d\n',fit.ChiSquare, fit.NDF)];
%     message = [message sprintf('   \\mu = %2.2f ps \\pm %2.2f ps\n',1000*mean_value,1000*fit.g1_mean_error)];
%     message = [message sprintf('   \\sigma_{1} = %2.2f ps \\pm %2.2f ps\n',1000*fit.sigma_core,1000*fit.sigma_core_error)];
%     message = [message sprintf('   \\sigma_{2} = %2.2f ps \\pm %2.2f ps\n',1000*fit.sigma_tail,1000*fit.sigma_tail_error)];
%     message = [message sprintf('   \\sigma_{comb} = %2.2f ps \\pm %2.2f ps \n',1000*fit.sigma_all,1000*fit.sigma_all_error)];
%     message = [message sprintf('   RMS_{tot} = %2.2f ps \\pm %2.2f ps', 1000*fit.StdDev, 1000*fit.StdDev/sqrt(2*fit.Entries))];

    xlabel('Time difference: PICOSEC vs reference, ns');
    ylabel('Events');
    legend('RAW hist','Gaussian 1', 'Gaussian 2','Gaussian comb');
    grid
    title_str = sprintf('Time resolution');
    title(title_str)
    xlim([-0.2 0.2])
    y_pos=get(gca,'ylim');
    x_pos=get(gca,'xlim');
    text(x_pos(1),0.8*y_pos(2),message)
    
    % Define the folder to save the PNG files
    save_folder = 'ResultsMATLAB';
    
    % Extract the run number from the file name
    [~, filename, ~] = fileparts(file_list_ROOT{file_index});
    run_number = regexp(filename, 'Run(\d+)_', 'tokens');
    run_number = run_number{1}{1};
    
    % Save the figure with the desired naming convention
    saveas(gcf, fullfile(save_folder, ['Run' run_number '_TimeResolution_2Gaussfit.png']));

end

