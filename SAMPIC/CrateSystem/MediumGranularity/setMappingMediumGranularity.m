%% multipad mapping

%format:
% MultipadPad SAMPICch
%format
%mmPadNumber SAMPIC channel
mappingArray = [];

% run July 25 run206
mappingTemp.mapping = [...
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
mappingTemp.run=206;

mappingArray=[mappingArray;mappingTemp];


% run July 25 run276
mappingTemp.mapping = [...
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
mappingTemp.run=276;

mappingArray=[mappingArray;mappingTemp];


% run July 25 run218
mappingTemp.mapping = [...
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

mappingTemp.run=218;
mappingArray=[mappingArray;mappingTemp];


runIDNumeric = str2num(run.id);

for mapPos = 1:length(mappingArray)
    mapEntry = mappingArray(mapPos);
    if mapEntry.run == runIDNumeric
        mapping = mapEntry.mapping;
        break;
    end
end

