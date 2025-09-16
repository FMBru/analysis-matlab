            for matchPos = 1:length(matchedEvents)
                    event = matchedEvents(matchPos);
                                                            str_disp=sprintf('Matched event %d', matchPos);
               % disp(str_disp);

                    matchedIdx = find(event.matchedChannels==triggerSignalsChannel);
                    if length(matchedIdx)==1
                        %found trigger match signal
                                        str_disp=sprintf('Found trig channel', i);
                disp(str_disp);

                        %go through all matched channels and analyse for the
                        %ones under test
                        
                        for dutChannelPos = 1:length(event.matchedChannels)
                            matchedIdxToAnalyse = find(channelToAnalyse==event.matchedChannels(dutChannelPos)); %check if this should be analysed
                            if length(matchedIdxToAnalyse)==1
                                %should analyse
                                                            str_disp=sprintf('Channel to analyse %d', event.matchedChannels(dutChannelPos));
                disp(str_disp);
                                %found match in DUT channel to analyse
                                eval(['dut = event.dut' num2str(event.matchedChannels(dutChannelPos)) ';']);
                                ref = event.refSignal;
                                %check match to SRS event counter
                                timeDiffTracker = scaledEventCounterTimes-ref.cell0Time;
                                [trackerMinValue,trackerMinIdx] = min(abs(timeDiffTracker));
                                trackerEntryMatched = eventCounterEntries(trackerMinIdx,:);
                                if size(trackerEntryMatched,1) %unambiguous match
                                    %found a single match in reference for this dut event -> unambiguois match
                                    if trackerMinValue<50000
                                        event.eventID = trackerEntryMatched(1);
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
                                        
                                        % process MM Picosec first to see if signal is valid
                                        if(minuit==1)
                                            %MM_temp = process_signal_minuit(t_vec_mm-etime,dutWaveform(:,2),opts_MM)
                                            MM_temp = process_signal_sampic(t_vec_mm-etime,dutWaveform(:,2),opts_MM,1);
                                        else
                                            MM_temp = process_signal(t_vec_mm-etime,dutWaveform(:,2),opts_MM);
                                        end
                                        
                                        % if signal is valid process MCP
                                        if(MM_temp.fail==0)
                                            if(minuit==1)
                                                %MCP_temp = process_signal_minuit(t_vec_mcp-etime,refWaveform(:,2),opts_MCP);
                                                MCP_temp = process_signal_sampic(t_vec_mcp-etime,refWaveform(:,2),opts_MCP,1);
                                            else
                                                MCP_temp = process_signal(t_vec_mcp-etime,refWaveform(:,2),opts_MCP);
                                            end
                                            
                                            % if MCP is valid store data to structure array
                                            if(MCP_temp.fail==0)
                                                % process tracker ID channel
                                                
                                                MM_temp.event_id = event_id;
                                                MCP_temp.event_id = event_id;
                                                
                                                
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
                                                
                                                % store valid data into array of structures
                                                MM_data(eventPosAnalysis)= MM_temp;
                                                MCP_data(eventPosAnalysis)= MCP_temp;
                                                time_diff(eventPosAnalysis) = MM_data(eventPosAnalysis).cfd.time-MCP_data(eventPosAnalysis).cfd.time;
                                                time_diff_sigmoid(eventPosAnalysis) = MM_data(eventPosAnalysis).sigmoid.timepoint-MCP_data(eventPosAnalysis).sigmoid.timepoint;
                                                MCP_maxy(eventPosAnalysis) = MCP_data(eventPosAnalysis).sig.max.y;
                                                MM_maxy(eventPosAnalysis) = MM_data(eventPosAnalysis).sig.max.y;
                                                trackerX(eventPosAnalysis) = MM_temp.x;
                                                trackerY(eventPosAnalysis) = MM_temp.y;
                                                eventIDArray(eventPosAnalysis) = MM_temp.event_id;
                                                dutChannelArray(eventPosAnalysis) = event.matchedChannels(dutChannelPos);
                                                
                                                eventPosAnalysis=eventPosAnalysis+1;
                                            end
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
                        
                       % matchedIdx = find(event.matchedChannels==channelToAnalyse);

                    end
                end