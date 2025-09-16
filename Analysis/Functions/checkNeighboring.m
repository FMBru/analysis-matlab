function [diffCol,diffRow] = isNeighboring(mainPad,testPad)
%mainPad = 25;
%testPad = 1;

colMain = mod(mainPad, 10);
rowMain = (mainPad - colMain) / 10; 

colTest = mod(testPad, 10);
rowTest = (testPad - colTest) / 10;

diffCol = abs(colTest-colMain);
diffRow = abs(rowTest-rowMain);

return;