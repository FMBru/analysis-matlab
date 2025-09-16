function isNeighbors = isNeighboring(mainPad,testPad)
%function [diffCol,diffRow] = isNeighboring(mainPad,testPad)
%mainPad = 25;
%testPad = 1;

colMain = mod(mainPad, 10);
rowMain = (mainPad - colMain) / 10; 

colTest = mod(testPad, 10);
rowTest = (testPad - colTest) / 10;

diffCol = abs(colTest-colMain);
diffRow = abs(rowTest-rowMain);

isNeighbors=1;

if diffCol>1
    isNeighbors=0;
end

if diffRow>1
    diffRow=0;
end

return;