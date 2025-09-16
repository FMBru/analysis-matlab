function [padMeanAmp] = moveAroundThePad(table, geom, padNum, store_folder, shouldSave, xlimits, ylimits)
%MOVEAROUNDTHEPAD Summary of this function goes here

figure
        
t = tiledlayout(3,3,'TileSpacing','none','Padding','none');
ax = gobjects(3,3);

for row = 1:3
    for col = 1:3
        x_centroid = geom.x_centroid - geom.x_window_length/3 + (col-1) * geom.x_window_length/3;
        y_centroid = geom.y_centroid - geom.y_window_length/3 + (row-1) * geom.x_window_length/3;
        
        x_first = x_centroid - geom.x_window_length/6;
        x_last = x_centroid + geom.x_window_length/6;
        %x_step = step;
        
        y_first = y_centroid - geom.y_window_length/6;
        y_last = y_centroid + geom.y_window_length/6;
        %y_step = step;
                
        geoCut = table.X > x_first & table.X < x_last & table.Y > y_first & table.Y < y_last;
        table_subset = table(geoCut, :);
        
        ax(row,col) = nexttile;
       

        message = correlationPlots(table_subset.e_peak_MM, table_subset.riseTime, 'e-peak-charge', 'pC', 'SAT', 'ns', store_folder, padNum, false, true);
        xlim(xlimits)
        ylim(ylimits)
        y_pos=get(gca,'ylim');
        x_pos=get(gca,'xlim');
        text(x_pos(1)+0.6*(x_pos(2)-x_pos(1)),y_pos(1)+0.7*(y_pos(2)-y_pos(1)),message, 'FontSize', 7);

        grid on;
        if row ~= 3
            ax(row,col).XTickLabel = [];
        end
        if col ~= 1
            ax(row,col).YTickLabel = [];
        end
        
    end
end

for col = 1:3
    linkaxes(ax(:,col), 'x');  
end

for row = 1:3
    linkaxes(ax(row,:), 'y');  
end

xlabel(t, 'e-peak charge, pC');
ylabel(t, 'Rise time, ns');

saveas(gcf,[store_folder '\riseTime_pad' num2str(padNum) '.png']);
end

