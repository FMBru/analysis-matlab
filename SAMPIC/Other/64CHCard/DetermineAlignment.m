nominalPad1X = 100;
nominalPad1Y = 100;
residuals = [];


for pos=1:length (padResults)
   padEntry = padResults(pos);
   
   padID = padEntry.padID
   
   row = floor(padID/10)+1;
   col = mod(padID,10);
   
   pad1OffsetXPads = (col-1)
   pad1OffsetYPads = (row-1)
   
   pad1OffsetXNominal = pad1OffsetXPads*10
   pad1OffsetYNominal = pad1OffsetYPads*10
   
   nominalPositionX = nominalPad1X-pad1OffsetXNominal;
   nominalPositionY = nominalPad1Y-pad1OffsetYNominal;
   
   if padID>0 &&  padEntry.xc~=0 && padEntry.xc~=100 && padEntry.xc~=25 && padEntry.yc~=0 && padEntry.yc~=100 && padEntry.yc~=25 
       xResidual = padEntry.xc-nominalPositionX;
       yResidual = padEntry.yc-nominalPositionY;
       residualEntry = [xResidual yResidual];
       residuals = [residuals;residualEntry];
       
       padEntry
       
         %  pause(5)

   end
    
    
end

close all
hist(residuals(:,1));
pause(2);
close all
hist(residuals(:,2));
pause(2);