%% multipad mapping

%format:
% MultipadPad SAMPICch

%please copy before changing

%format
%mmPadNumber SAMPIC channel

mappingArray = [];

% run July 25 run206
mapTemp.mapping = [...
    1 1;
    2 20;
    3 3;
    4 4;
    5 5;
    6 6;
    7 7;
    8 8;
    9 9; 
    10 10;
    11 11;
    12 12 ;
    13 13;
    14 14;
    15 15;
    16 16; 
    17 21;
    18 18;
    19 19;
];
mapTemp.runID = 206;
mappingArray = [mappingArray;mapTemp];

% run July 25 run276
mapTemp.mapping = [...
    1 1;
    2 2;
    3 3;
    4 4;
    5 5;
    6 6;
    7 22;
    8 8;
    9 9; 
    10 10;
    11 11;
    12 12 ;
    13 13;
    14 14;
    15 15;
    16 16; 
    17 20;
    18 18;
    19 19;
];
mapTemp.runID = 276;
mappingArray = [mappingArray;mapTemp];

% run July 25 run218
mapTemp.mapping = [...
    1 1;
    2 20;
    3 3;
    4 4;
    5 5;
    6 6;
    7 7;
    8 8;
    9 9; 
    10 22;
    11 11;
    12 12 ;
    13 13;
    14 14;
    15 15;
    16 16; 
    17 21;
    18 18;
    19 19;
];

mapTemp.runID = 218;
mappingArray = [mappingArray;mapTemp];

%select mapping based on runID

for mappingPos = 1:length(mappingArray)
    mapEntry = mappingArray(mappingPos);
    if str2num(run.id) == mapEntry.runID
        mapping = mapEntry.mapping;
        break;
    end
end
