function table = cleanDouble(table)
%CLEANDOUBLE clean double eventID in sorted table
table(ismember(table.eventID,table.eventID(diff(table.eventID) == 0)), :) = [];
end

