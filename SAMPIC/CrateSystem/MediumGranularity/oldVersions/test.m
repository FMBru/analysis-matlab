for i = 1:19
    filtered_data = padAmplitudes{i};
    % Export data to .txt for later analysis
    outputDir = [store_folder '/AmpExtractedData'];
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end
    fullPath = [outputDir '/AmpPad' num2str(i) '.txt'];
    fileID = fopen(fullPath, 'w');
    fprintf(fileID, '%f\n', filtered_data);
    fclose(fileID);
    
    if ~isempty(filtered_data)
        pause(1);
        figure;
        h = histogram(filtered_data, 50);
        
        % Filter out empty bins for fitting
        bin_centers = h.BinEdges(1:end-1) + h.BinWidth / 2;
        non_zero_bins = h.Values > 0;
        
        
            % Fit the curve using only non-zero bins
            [fitted_curve, gof] = fit(bin_centers(non_zero_bins)', h.Values(non_zero_bins)', fitfun, options);
            
            % Store the fitted parameters
            P = [P; fitted_curve.N, fitted_curve.theta, fitted_curve.nBar];
            
            % Plot the fitted curve
            if fitted_curve.theta > 0
                hold on;
                plot(bin_centers(non_zero_bins), fitted_curve(bin_centers(non_zero_bins)), 'LineWidth', 3, 'Color', [0, 0, 1]);
                
                % Calculate chi-squared and degrees of freedom
                ch2 = sum(((h.Values(non_zero_bins) - fitted_curve(bin_centers(non_zero_bins))').^2) ./ fitted_curve(bin_centers(non_zero_bins))');
                dof = sum(non_zero_bins) - 3;
                nch2 = ch2 / dof;
                np = 1 - chi2cdf(ch2, dof); % P(chi^2 > ch2)
                
                rmssAmplitudes = [rmssAmplitudes, rms(filtered_data)];
                musAmplitudes = [musAmplitudes, fitted_curve.nBar];
        
            end
            
            % Display results
            message = sprintf('amplitude = %2.1f mV', 1000*fitted_curve.nBar);
            y_pos = get(gca, 'ylim');
            x_pos = get(gca, 'xlim');
            text(x_pos(1), 0.75 * y_pos(2), message);
            xlabel('Amp (V)');
            ylabel('Events');
            title_str = sprintf('%s \n Amp hist  - Pad %d', runTitleString, i);
            title(title_str);
            grid on;
            
            % Save the figure
            saveas(gcf, [store_folderPads '\AmpHist-Run' run.id '_pad' num2str(i) '.png']);
            close all;
        end
    else
        musAmplitudes = [musAmplitudes, 0];
        sigmasAmplitudes = [sigmasAmplitudes, 0];
        rmssAmplitudes = [rmssAmplitudes, 0];
    end
end

VisualiseMediumGranularity(channelsEnabled, musAmplitudes, 'Map Multipad Picosec Mean Amplitude', [store_folder '\Run' run.id '_AmpMAP.png'], 'Amp (V)', '%0.5f %', 0.03, 0.06);

