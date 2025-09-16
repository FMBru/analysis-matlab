function [out] = diagPlots(table,pad_geom, quantity, padNumber, store_folder, shouldSave)
%DIAGPLOTS 

nSteps = 7;

step = pad_geom.pad_window_length / nSteps;
x_points = (pad_geom.x_centroid - pad_geom.pad_window_length/2 + step/2):step:(pad_geom.x_centroid + pad_geom.pad_window_length/2 - step/2);
y_points = (pad_geom.y_centroid - pad_geom.pad_window_length/2 + step/2):step:(pad_geom.y_centroid + pad_geom.pad_window_length/2 - step/2);
y_points_inverted = y_points(end:-1:1);

cut_radius = step;
bottom_left_vec = zeros(1,nSteps);
bottom_left_vec_err = zeros(1,nSteps);
top_left_vec = zeros(1,nSteps);
top_left_vec_err = zeros(1,nSteps);

for j=1:nSteps

    bottom_left_circ = (table.X - x_points(j)).^2 + (table.Y - y_points(j)).^2 < cut_radius^2;
    subset = table(bottom_left_circ, :);
    bottom_left_vec(j) = mean(subset.(quantity));
    bottom_left_vec_err(j) = std(subset.(quantity)) / sqrt(height(subset));

    top_left_circ = (table.X - x_points(j)).^2 + (table.Y - y_points_inverted(j)).^2 < cut_radius^2;
    subset = table(top_left_circ, :);
    subset = subset(~isnan(subset.(quantity)),:);
    top_left_vec(j) = mean(subset.(quantity));
    top_left_vec_err(j) = std(subset.(quantity)) / sqrt(height(subset));
    
end

figure;
errorbar(x_points, top_left_vec, top_left_vec_err, top_left_vec_err, 'o', 'LineWidth',1);
hold on
errorbar(x_points, bottom_left_vec, bottom_left_vec_err, bottom_left_vec_err, 'o', 'LineWidth',1);
xlabel('Distance from pad center, mm')
ylabel(quantity)
grid on
title([quantity ' along the diagonals of pad ' num2str(padNumber)]);
legend('Top left diag.', 'Bottom left diag.');
if shouldSave
    saveas(gcf,[store_folder '\pad' num2str(padNumber) '\' quantity 'diagonalDistribution_pad' num2str(padNumber) '.png']);
end
close

