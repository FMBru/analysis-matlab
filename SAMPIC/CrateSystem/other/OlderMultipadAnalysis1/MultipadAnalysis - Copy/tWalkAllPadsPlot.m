% Define the directory containing the tWalkData files and correction files
close all
clear all
runIDsToAnalyse = ["232"];
run.id = convertStringsToChars(runIDsToAnalyse(1));
store_folder = ['\\eosproject-smb\eos\project\p\picosec\testbeam\2024_June_h4\Results\Run' run.id '-SAMPIC'];
dataDir = [store_folder '/tWalkData/'];
folder = ['\\eosproject-smb\eos\project\p\picosec\testbeam\2024_June_h4\Results\Run232-SAMPIC\Pads\'];

% Get a list of all tWalkData files in the directory
filePattern = fullfile(dataDir, 'tWalkData_Pad*.txt');
dataFiles = dir(filePattern);

% Check if any files are found
if isempty(dataFiles)
    error('No tWalkData files found in directory: %s', dataDir);
end

% Manually define a large set of distinct colors
colors = [
    0 0.4470 0.7410; % Blue
    0.8500 0.3250 0.0980; % Red
    0.9290 0.6940 0.1250; % Yellow
    0.4940 0.1840 0.5560; % Purple
    0.4660 0.6740 0.1880; % Green
    0.3010 0.7450 0.9330; % Cyan
    0.6350 0.0780 0.1840; % Dark Red
    0.8 0.4 0.6; % Custom Pink
    0.2 0.6 0.8; % Custom Light Blue
    0.3 0.3 0.3; % Grey
    0.9 0.7 0.3; % Custom Gold
    0.7 0.7 0.2; % Olive
    0.4 0.6 0.4; % Forest Green
    0.9 0.5 0.5; % Coral
    0.5 0.5 0.7; % Slate Blue
    0.6 0.6 0.2; % Mustard
];
numColors = length(colors);

% Pads to skip
badPads = [9, 16, 17];
latency = [28.4, 28.1, 28.4, 28.3, 28.1, 26.9, 28.4, 27.6, 26.1, 28.1, 28, 27, 26.8, 25.6, 26.1, 25.3, 0, 27.6, 26.9];
satOffset = zeros(19,1);

% Create a figure with subplots for the marginal histograms
figure;

% Define the relative sizes of the main plot and histograms
mainPlotPos = [0.15 0.15 0.6 0.6]; % Main plot position
histHeight = 0.2; % Height of the histogram plots
histWidth = 0.2; % Width of the histogram plots

% Main plot axes
mainPlot = axes('Position', mainPlotPos);
xlabel('Electron peak amplitude (V)');
ylabel('SAT, ps');
grid on;
hold on;

% Initialize variables for stacking histograms
numBins = 50;  % Define the number of bins
xBinEdges = linspace(0, 0.4, numBins + 1);
yBinEdges = linspace(-1, 1, numBins + 1);
xHistCounts = zeros(numBins, 1);
yHistCounts = zeros(numBins, 1);

% Loop over each file and plot the data
for k = 1:length(dataFiles)
    % Extract the pad number from the filename
    baseFileName = dataFiles(k).name;
    padNumber = sscanf(baseFileName, 'tWalkData_Pad%d.txt');
    
    % Skip the bad pads
    if ismember(padNumber, badPads)
        continue;
    end

    % Read the data from the file
    fullFileName = fullfile(dataFiles(k).folder, baseFileName);
    data = load(fullFileName);
    
    % Construct the correction file name
    correctionFileName = sprintf('twCorrection_Pad%d.txt', padNumber);
    correctionFilePath = fullfile(folder, correctionFileName);
    
    % Check if the correction file exists
    if exist(correctionFilePath, 'file') == 2
        % Read the value to subtract from the correction file
        correctionData = load(correctionFilePath);
        value_to_subtract = correctionData(1, 1); % Extract the first value
        disp(['Value to subtract for Pad ', num2str(padNumber), ': ', num2str(value_to_subtract)]);
    else
        error('Correction file does not exist: %s', correctionFilePath);
    end
    
    % The data is assumed to be in the format:
    % e_peak, mean_sat, e_peak_err, e_peak_err_n, e_peak_err_p
    e_peak = data(:, 1);

    % Determine the index where the last 30% begins
    %startIndex = ceil(length(data(:, 2)) * 0.7);
    %correction = median(data(startIndex:end, 2));

    mean_sat = data(:, 2) - value_to_subtract;%correction; % Adjust mean_sat
    e_peak_err = data(:, 3);

    % Plot the data with error bars and assign a unique color
    colorIndex = mod(k-1, numColors) + 1; % Cycle through colors if needed
    errorbar(mainPlot, e_peak, mean_sat, e_peak_err/2, 'o', ...
             'MarkerSize', 2, 'MarkerFaceColor', colors(colorIndex,:), ...
             'CapSize', 0, 'LineWidth', 1, 'Color', colors(colorIndex,:), ...
             'DisplayName', ['Pad ' num2str(padNumber)]);

   %plot(mainPlot,e_peak,1000*(twalk_fn_minuit(correctionData(1,:),e_peak')-correctionData(1,1)),'LineWidth',1.3, 'Color', colors(colorIndex,:), 'DisplayName', ['Pad ' num2str(padNumber)]);


    % Store the histogram counts for stacking
    xHistCounts = xHistCounts + histcounts(e_peak, xBinEdges)';
    yHistCounts = yHistCounts + histcounts(mean_sat, yBinEdges)';
end

% Plot stacked histograms with colors
% X-axis marginal histogram (on top)
xHistPos = [mainPlotPos(1) mainPlotPos(2) + mainPlotPos(4) + 0.05 mainPlotPos(3) histHeight];
xHistPlot = axes('Position', xHistPos);
bar(xBinEdges(1:end-1), xHistCounts, 'FaceColor', 'flat', 'EdgeColor', 'none');
colormap(xHistPlot, colors);
set(xHistPlot, 'XTick', [], 'YColor', 'k'); % Hide x-axis ticks, show y-axis ticks
ylabel('Number of Events');
axis tight;

% Y-axis marginal histogram (on the right)
yHistPos = [mainPlotPos(1) + mainPlotPos(3) + 0.05 mainPlotPos(2) histWidth mainPlotPos(4)];
yHistPlot = axes('Position', yHistPos);
barh(yBinEdges(1:end-1), yHistCounts, 'FaceColor', 'flat', 'EdgeColor', 'none');
colormap(yHistPlot, colors);
set(yHistPlot, 'YTick', [], 'XColor', 'k'); % Hide y-axis ticks, show x-axis ticks
xlabel('Number of Events');
axis tight;

% Link the main plot with the marginal histograms
linkaxes([mainPlot xHistPlot], 'x');
linkaxes([mainPlot yHistPlot], 'y');

% Set limits for the axes
xlim(mainPlot, [0 0.4]); % Example x-axis limit, adjust as needed
ylim(mainPlot, [-0.5 1]); % Example y-axis limit, adjust as needed

% Show the legend only on the main plot
legend(mainPlot, 'show');

% Move the plot window to a specified location on the screen
movegui(gcf, 'center');
