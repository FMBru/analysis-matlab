function T_clean = removeAllNaNRows(T, keepCols)
%REMOVEALLNANROWS Removes events in which all entries (except in keepCols) are NaN

    allCols = T.Properties.VariableNames;
    checkCols = setdiff(allCols, keepCols);

    allNaN = all(isnan(T{:, checkCols}), 2);

    T_clean = T(~allNaN, :);
end