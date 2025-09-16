                        matchedEvents = [];

            for refPos = 1:refSignalsStep:length(signalsRefSegment)
                refSignal = signalsRefSegment(refPos);
                
                % Go through all channels and check if there is a match
                matchedChannels = [];
                
                for chDUTPos = 1:length(channelsEnabled)
                    if  channelsEnabled(chDUTPos) ~= referenceDetectorChannel
                        shouldAnalyseChannel = find(channelsEnabled(chDUTPos)==channelToAnalyse);
                        if length(shouldAnalyseChannel)==1
                            %is dut channel
                            latency = chLatencyArray(chDUTPos);
                            dutChannelMask = chIDsArraySegment==channelsEnabled(chDUTPos);
                            chTimesArrayDUT = chTimesArraySegment(dutChannelMask);
                            chIndexArrayDUT = chIndexArraySegment(dutChannelMask);
                            
                            minTime = latency-matchingDUTREFWindowWidth;
                            maxTime = latency+matchingDUTREFWindowWidth;
                            timeDiffArray = abs(chTimesArrayDUT-refSignal.cell0Time);
                            
                            %minTimeDiff = min(timeDiffArray);
                            %minTimeDiffs = [minTimeDiffs;minTimeDiff];
                            
                            matchMask = timeDiffArray>minTime & timeDiffArray<maxTime;
                            
                            %chIndexMatched = chIndexArrayDUT(matchMask); %wrong?
                            chIndexMatched = timeDiffArray(matchMask);
                            
                            if length(chIndexMatched)==1 %unambiguous hit
                                %accept as match
                                chIndex = chIndexArrayDUT(matchMask);
                                
                                %time between signals
                                signalTemp = signalsSegment(chIndex);
                                
                                correctTime = timeDiffArray(matchMask)
                                timeDiffTemp = (signalTemp.cell0Time-refSignal.cell0Time)
                                
                         if timeDiffTemp>100
                            pause(5);
                        end                                

                                timesTestArray = [timesTestArray;[correctTime timeDiffTemp]];
                                
                                
                                %dutTemp = signalsSegment(chIndexMatched(1));
                                dutTemp = signalTemp;
                                if length(matchedChannels)==0
                                    matchedChannels = [channelsEnabled(chDUTPos)];
                                    
                                else
                                    matchedChannels(length(matchedChannels)+1) = channelsEnabled(chDUTPos);
                                    
                                end
                                
                            else
                                dutTemp = [];
                            end
                            
                            eval(['event.dut' num2str(channelsEnabled(chDUTPos)) ' = dutTemp;']);
                        end
                    end
                    
                end
                event.matchedChannels = matchedChannels;
                event.numberMatched = length(matchedChannels);
                event.refSignal = refSignal;
                nextIndex = length(matchedEvents)+1;
                
%                 if event.numberMatched > 0
%                     event      
%                     pause(1);
%                 end

                if event.numberMatched>1
                    if nextIndex==1
                        matchedEvents = [event];
                    else
                        
                       
                       
                        matchedEvents(nextIndex) = event;
                        

                    end
                end
            end