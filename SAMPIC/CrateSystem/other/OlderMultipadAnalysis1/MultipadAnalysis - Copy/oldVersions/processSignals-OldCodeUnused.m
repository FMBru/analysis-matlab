if length(signals)>5000
                %% match and process signals
                str_disp=sprintf('Matching events', i);
                disp(str_disp);
                matchedEvents = [];

                minSignalsTime = min(chTimesArray);
                maxSignalsTime = max(chTimesArray);

                chTimesArrayRef = chTimesRefArray;
                chIndexArrayRef = chIndexRefArray;
                refSignalsSubsetMask = chTimesArrayRef>minSignalsTime & chTimesArrayRef<maxSignalsTime;

                signalsRefSubset = signalsRef(refSignalsSubsetMask);
