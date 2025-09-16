% Define the path for the MAT file
storeMatfileFolder = 'C:\Users\GDD\Documents\Picosec\Run292-SAMPIC\variables';
processingBatchCounter = 13;
storeMatfilePath = [storeMatfileFolder '\processedBatch_' int2str(processingBatchCounter) '.mat'];

% Load the data from the MAT file
loadedData = load(storeMatfilePath);

% Extract the relevant data fields into variables
MM_data = loadedData.MM_data;
MCP_data = loadedData.MCP_data;
time_diff = loadedData.time_diff;
time_diff_sigmoid = loadedData.time_diff_sigmoid;
MCP_maxy = loadedData.MCP_maxy;
MM_maxy = loadedData.MM_maxy;
trackerX = loadedData.trackerX;
trackerY = loadedData.trackerY;
eventIDArray = loadedData.eventIDArray;
dutChannelArray = loadedData.dutChannelArray;
currentChannel = loadedData.currentChannel;
MM_bgAvg = loadedData.MM_bgAvg;
MM_bgRMS = loadedData.MM_bgRMS;

% Create a table to hold all the data
dataTable = table(MM_data, MCP_data, time_diff, time_diff_sigmoid, MCP_maxy, MM_maxy, ...
                  trackerX, trackerY, eventIDArray, dutChannelArray, currentChannel, ...
                  MM_bgAvg, MM_bgRMS);
% dataTable = table(MM_data.fail, MM_data.cfd,MM_data.sig, MM_data.sigmoid, MM_data.event_id, MM_data.waveform, MM_data.x, MM_data.y, MM_data.sampic);             
% Define the path for the CSV file
storeCsvFilePath = [storeMatfileFolder 'MM_data_test.csv'];

% Write the table to a CSV file
writetable(dataTable, storeCsvFilePath, 'WriteVariableNames', false);

% Display a message indicating the CSV file was saved
disp(['Data saved to ' storeCsvFilePath]);

% Plotting example data
figure;
plot(MM_data(1).waveform(:, 2));
