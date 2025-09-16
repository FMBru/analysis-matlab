clear all

addpath 'C:\Users\GDD\Documents\MATLAB\Picosec\SAMPIC\CrateSystem\other'

%%define variables necessary for analyses

dutChannelArray = [];
time_diff_sigmoid = [];
MCP_maxy = [];
MM_maxy = [];
trackerX = [];
trackerY = [];
eventIDArray = [];
MM_epeak = [];
MCP_epeak = [];
MM_riseTime = [];
MM_bgAvg = [];
MM_bgRMS = [];
MM_sampicTOT = [];
MM_sampicAmp = [];

%%enter run IDs to analyse together
runIDsToAnalyse = ["417"];

for runPos = 1:length(runIDsToAnalyse)

run.id = convertStringsToChars(runIDsToAnalyse(runPos));

run.name = run.id;
run.oscilloscope = 'SAMPIC';


%variablesFolder = ['\\eosproject-smb\eos\project\p\picosec\testbeam\2022_October_h4\Results\Run' run.id '-SAMPIC\variables'];
variablesFolder = ['C:\Users\GDD\Documents\Picosec\Run' run.id '-SAMPIC\variables'];
str_disp=sprintf('Loading variables for run %s', run.id);
disp(str_disp);
variablesFilesList = dir(variablesFolder);
numberFiles = length(variablesFilesList);
str_disp=sprintf('Found %d files to load', numberFiles);
disp(str_disp);

for i=1:length(variablesFilesList)
    fileInfo = variablesFilesList(i);
    
    str_disp=sprintf('Loading file %d / %d', i, numberFiles);
    disp(str_disp);
    
    if contains(fileInfo.name,'.mat')
        filePath = [variablesFolder '\' fileInfo.name];
        segment = load(filePath);
        numberEventsLoadedSegment = length(segment.time_diff_sigmoid);
        numberEventsLoadedSegment = length(segment.MM_data);
        
        dutChannelArray = [dutChannelArray segment.dutChannelArray];
        time_diff_sigmoid = [time_diff_sigmoid segment.time_diff_sigmoid];
        MCP_maxy = [MCP_maxy segment.MCP_maxy];
        MM_maxy = [MM_maxy segment.MM_maxy];
        trackerX = [trackerX segment.trackerX];
        trackerY = [trackerY segment.trackerY];
        eventIDArray = [eventIDArray segment.eventIDArray];
        MM_bgAvg = [MM_bgAvg segment.MM_bgAvg];
        MM_bgRMS = [MM_bgRMS segment.MM_bgRMS];

        for pos = 1:numberEventsLoadedSegment
            %step through loaded file and assemble arrays
            MM_epeakTemp(pos) = segment.MM_data(pos).sig.charge.e_peak;
            MCP_epeakTemp(pos)= segment.MCP_data(pos).sig.charge.e_peak;% extract electrom peak charge
            MM_riseTimeTemp(pos)= segment.MM_data(pos).sigmoid.timepoint90- segment.MM_data(pos).sigmoid.timepoint10;% extract electrom peak charge
            MM_sampicTOTTemp(pos)= segment.MM_data(pos).sampic.tot;
            MM_sampicAmpTemp(pos)= -segment.MM_data(pos).sampic.amp;
        end
        MM_epeak = [MM_epeak MM_epeakTemp];
        MCP_epeak = [MCP_epeak MCP_epeakTemp];
        MM_riseTime = [MM_riseTime MM_riseTimeTemp];
        MM_sampicTOT = [MM_sampicTOT MM_sampicTOTTemp];
        MM_sampicAmp = [MM_sampicAmp MM_sampicAmpTemp];
       
        clear MM_epeakTemp MCP_epeakTemp MM_riseTimeTemp MM_sampicTOTTemp MM_sampicAmpTemp
        clear segment
    end
end

    str_disp=sprintf('Loaded %d events', length(dutChannelArray));
    disp(str_disp);
end

AnalyseRunSAMPIC_TOT
