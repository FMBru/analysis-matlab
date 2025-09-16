rmsResults = [];

for resultPos=1:length(padResults)
        if padResults(resultPos).padID>0 && padResults(resultPos).numberEntriesSampling>900 && padResults(resultPos).rmsSampling>0 && padResults(resultPos).numberEntriesSampling>20

    rmsResults = [rmsResults;padResults(resultPos).rmsSampling];
        end
end
figure
h = histogram(rmsResults,15);
xlabel('Time resolution RMS, ps');
ylabel('Number pads');
xlim([10 40]);
ylim([0 12]);
grid on