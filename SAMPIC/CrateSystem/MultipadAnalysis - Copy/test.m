
            %have matched events - go through and process events
            str_disp=sprintf('Processing events', i);
            disp(str_disp);
            eventPosProcessing = 1;
            
            for matchPos = 1:matchedEventsStep:length(matchedEvents)
                
                if mod(matchPos,5000)==0
                    toc
                    tic
                    str_disp=sprintf('Processing event No. %d', matchPos);
                    disp(str_disp);
                end
                
                event = matchedEvents(matchPos);
                for dutChannelPos = 1:length(event.matchedChannels)
                    matchedIdxToAnalyse = find(channelToAnalyse==event.matchedChannels(dutChannelPos)); %check if this should be analysed
                    if length(matchedIdxToAnalyse)==1
                        %should analyse
                        eval(['dut = event.dut' num2str(event.matchedChannels(dutChannelPos)) ';']);
                        ref = event.refSignal;
                        %check match to SRS event counter
                        timeDiffTracker = scaledEventCounterTimes-ref.cell0Time;
                        [trackerMinValue,trackerMinIdx] = min(abs(timeDiffTracker));
                        trackerEntryMatched = eventCounterEntries(trackerMinIdx,:);
                        if size(trackerEntryMatched,1) %unambiguous match
                            %found a single match in reference for this dut event -> unambiguois match
                            trackerMinValue;
                            if trackerMinValue<1e15
                                % Double defined to avoid printing
                                eID = trackerEntryMatched(1);
                                event.eventID = eID;
                                %matched to tracker hit, process event
                                dutWaveform = dut.waveform;
                                refWaveform = ref.waveform;
                                event_id = event.eventID;
                                
                                % generate virtual time vector (important to take care for trigger
                                % offset in LeCroy scope)
                                t_vec_mm=dutWaveform(:,1);
                                t_vec_mcp=refWaveform(:,1);
                                
                                % subtract the earliest time
                                etime=min([t_vec_mm; t_vec_mcp]);
                                run.savedSignals = run.savedSignals +1;
                                if run.savedSignals<2
                                    shouldSave = 1;
                                else
                                    shouldSave = 0;
                                end
                                baselineDUT = mean(dutWaveform(1:5,2));
                                baselineMCP = mean(refWaveform(1:5,2));
                                
                                % process MM Picosec first to see if signal is valid
                                MM_temp = process_signal_sampic(t_vec_mm-etime,dutWaveform(:,2),opts_MM,baselineDUT,run,shouldSave,storeFolderSignalDUT);
                                
                                % if signal is valid process MCP
                                if(MM_temp.fail==0)
                                    MCP_temp = process_signal_sampic(t_vec_mcp-etime,refWaveform(:,2),opts_MCP,baselineMCP,run,shouldSave,storeFolderSignalRef);
                                    
                                    % if MCP is valid store data to structure array
                                    if(MCP_temp.fail==0)
                                        % process tracker ID channel
                                        
                                        MM_temp.event_id = event_id;
                                        MCP_temp.event_id = event_id;
                                        
                                        MM_temp.waveform = [(t_vec_mm-etime) (dutWaveform(:,2))];
                                        MCP_temp.waveform = [(t_vec_mm-etime) (refWaveform(:,2))];
                                        
                                        event_id;
                                        xPos = 0;
                                        yPos = 0;
                                        if tracker.en
                                            %find corresponding tracker entry and save XY info
                                            trackIdx = find(tracker.data(:,1) == event_id);
                                            if trackIdx>0
                                                %is valid trackIdx -> get XY
                                                xPos = tracker.data(trackIdx,dutColArray(tracker.dutIndex,2));
                                                yPos = tracker.data(trackIdx,dutColArray(tracker.dutIndex,3));
                                                if size(xPos,1)>1
                                                    xPos = xPos(1);
                                                end
                                                if size(yPos,1)>1
                                                    yPos = yPos(1);
                                                end
                                            end
                                            
                                        end
                                        MM_temp.x = xPos;
                                        MM_temp.y = yPos;
                                        
                                        %save data directly extracted sampic file
                                        %in MM_temp (amp, TOT)
                                        
                                        MM_temp.sampic.tot = dut.TOTvalue;
                                        MM_temp.sampic.amp = dut.amplitude;
                                        
                                        MCP_temp.sampic.tot = ref.TOTvalue;
                                        MCP_temp.sampic.amp = ref.amplitude;
                                        
                                        % store valid data into array of structures
                                        MM_data(eventPosProcessing)= MM_temp;
                                        MCP_data(eventPosProcessing)= MCP_temp;
                                        time_diff(eventPosProcessing) = MM_data(eventPosProcessing).cfd.time-MCP_data(eventPosProcessing).cfd.time;
                                        time_diff_sigmoid(eventPosProcessing) = MM_data(eventPosProcessing).sigmoid.timepoint-MCP_data(eventPosProcessing).sigmoid.timepoint;
                                        MCP_maxy(eventPosProcessing) = MCP_data(eventPosProcessing).sig.max.y;
                                        MM_maxy(eventPosProcessing) = MM_data(eventPosProcessing).sig.max.y;
                                        MM_bgAvg(eventPosProcessing) = MM_data(eventPosProcessing).sig.blavg;
                                        MM_bgRMS(eventPosProcessing) = MM_data(eventPosProcessing).sig.blrms;
                                        trackerX(eventPosProcessing) = MM_temp.x;
                                        trackerY(eventPosProcessing) = MM_temp.y;
                                        eventIDArray(eventPosProcessing) = MM_temp.event_id;
                                        dutChannelArray(eventPosProcessing) = event.matchedChannels(dutChannelPos);
                                        currentChannel = event.matchedChannels(dutChannelPos);
                                        
                                        eventPosProcessing=eventPosProcessing+1;
                                        
                                        %                                                 %plot signals as examples
                                        %                                                 %of events
                                        %                                                 plot(t_vec_mm-etime,ch_mm.y(:,m));
                                        %                                                 hold on;
                                        %                                                 plot(t_vec_mcp-etime,ch_mcp.y(:,m));
                                        %                                                 grid on
                                        %                                                 xlim([200e-9 300e-9])
                                        %                                                 title('DUT / REF');
                                        %                                                 ylabel('Voltage, V');
                                        %                                                 xlabel('Time, ns');
                                        %                                                 legend('DUT','REF');
                                        %                                                 grid on;
                                    else
                                        disp('MCP processing failed');
                                        
                                    end
                                else
                                end
                                
                            else
                                
                                
                                event.eventID = 0;
                            end
                            %eventPos = eventPos+1;
                        else
                            %found multiple matches
                        end
                        
                    end
                end
                
                
                if eventPosProcessing>10000
                    %save batch
                    str_disp=sprintf('Processed batch, total events analysed: %i latestEventID: %i', eventPosProcessing,eventIDArray(eventPosProcessing-1));
                    disp(str_disp);
                    n_strPadded = sprintf( '%03d', processingBatchCounter ) ;                    
                    storeMatfilePath = [storeMatfileFolder '\processedBatch_' n_strPadded '.mat'];
                    m = matfile(storeMatfilePath,'Writable',true);
                    %
                    save(storeMatfilePath,'MM_data');
                    save(storeMatfilePath,'MCP_data','-append');
                    save(storeMatfilePath,'time_diff','-append');
                    save(storeMatfilePath,'time_diff_sigmoid','-append');
                    save(storeMatfilePath,'MCP_maxy','-append');
                    save(storeMatfilePath,'MM_maxy','-append');
                    save(storeMatfilePath,'trackerX','-append');
                    save(storeMatfilePath,'trackerY','-append');
                    save(storeMatfilePath,'eventIDArray','-append');
                    save(storeMatfilePath,'dutChannelArray','-append');
                    save(storeMatfilePath,'currentChannel','-append');
                    save(storeMatfilePath,'MM_bgAvg','-append');
                    save(storeMatfilePath,'MM_bgRMS','-append');
                    
                    clear MM_data MCP_data time_diff time_diff_sigmoid MCP_maxy MM_maxy trackerX trackerY eventIDArray dutChannelArray currentChannel MM_bgAvg MM_bgRMS;
                    eventPosProcessing = 1;
                    processingBatchCounter = processingBatchCounter+1;
                    str_disp=sprintf('Saved variables');
                    disp(str_disp);
                    
                end
                
                
            end
            
            
        
    