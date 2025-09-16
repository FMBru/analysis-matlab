clear all

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

%%enter run IDs to analyse together
runIDsToAnalyse = ["206"];

for runPos = 1:length(runIDsToAnalyse)

run.id = convertStringsToChars(runIDsToAnalyse(runPos));

run.name = run.id;
run.oscilloscope = 'SAMPIC';


%variablesFolder = ['\\eosproject-smb\eos\project\p\picosec\testbeam\2022_October_h4\Results\Run' run.id '-SAMPIC\variables'];
    variablesFolder = ['F:\Processing\Picosec\Run' run.id '-SAMPIC\variables'];
str_disp=sprintf('Loading variables for run %s', run.id);
disp(str_disp);
variablesFilesList = dir(variablesFolder);
numberFiles = length(variablesFilesList);
str_disp=sprintf('Found %d files to load', numberFiles);
disp(str_disp);

for i=1:length(variablesFilesList)
    fileInfo = variablesFilesList(i);
    
    
    if contains(fileInfo.name,'.mat') & contains(fileInfo.name,'processedBatch')
        
           str_disp=sprintf('Loading file %d / %d', i, numberFiles);
           disp(str_disp);

    
        filePath = [variablesFolder '\' fileInfo.name];
        segment = load(filePath);
        numberEventsLoadedSegment = length(segment.time_diff_sigmoid);
        
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
        end
        MM_epeak = [MM_epeak MM_epeakTemp];
        MCP_epeak = [MCP_epeak MCP_epeakTemp];
        MM_riseTime = [MM_riseTime MM_riseTimeTemp];
        clear MM_epeakTemp MCP_epeakTemp MM_riseTimeTemp
        
        clear segment
    end
end

    str_disp=sprintf('Loaded %d events', length(dutChannelArray));
    disp(str_disp);
end

AnalyseAllChannels
