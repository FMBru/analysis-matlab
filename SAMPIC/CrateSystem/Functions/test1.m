padIDsArray = zeros(10,10);
padMeanArray = zeros(10,10);

for row = 1:10
    for col=1:10
        padID = 10*(row-1)+col;
        padIDsArray (row,col) = padID;
    end
end