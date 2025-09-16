function [pad_x_center,pad_y_center, pad_mask] = findPadCenter(geom_frame, map, pad_window_length, tol)
%FINDPADCENTER' this function tries to find the pad center: first it looks
%at the sum of the amplitudes of the singlas to look for a naive peak and
%then makes a cut around this peak to fit the projections and find the mean
%(pad_mask at the end is a mask that is tol*pad_window_length)

xCenters = geom_frame.xCenters;
yCenters = geom_frame.yCenters;
Xgrid = geom_frame.Xgrid;
Ygrid = geom_frame.Ygrid;

projX = sum(map, 1, 'omitnan');   % sum over rows → 1 × numXBins
projY = sum(map, 2, 'omitnan');   % sum over columns → numYBins x 1

figure;
subplot(1,2,1);
bar(xCenters, projX);
xlabel('X'); ylabel('Sum amplitude'); title('Projection on X');

subplot(1,2,2);
bar(yCenters, projY);
xlabel('Y'); ylabel('Sum amplitude'); title('Projection on Y');

hold on

%setting the limits for pad center fit
[~, pad_center_naive_x_idx] = max(projX);
pad_center_naive_x = xCenters(pad_center_naive_x_idx);
[~, pad_center_naive_y_idx] = max(projY);
pad_center_naive_y = yCenters(pad_center_naive_y_idx);
pad_window_left = pad_center_naive_x - pad_window_length/2;
pad_window_right = pad_center_naive_x + pad_window_length/2;
pad_window_up = pad_center_naive_y + pad_window_length/2;
pad_window_bottom = pad_center_naive_y - pad_window_length/2;

subplot(1,2,1);
xline(pad_window_left, '--r', 'LineWidth', 2);
xline(pad_window_right, '--r', 'LineWidth', 2);
xline(pad_center_naive_x, '--g', 'LineWidth', 2);

subplot(1,2,2);
xline(pad_window_up, '--r', 'LineWidth', 2);
xline(pad_window_bottom, '--r', 'LineWidth', 2);
xline(pad_center_naive_y, '--g', 'LineWidth', 2);

% isolate the center of the pad and use a fit  to find the center

pad_mask = Xgrid <= pad_window_right & Xgrid >= pad_window_left & Ygrid <= pad_window_up & Ygrid >= pad_window_bottom;
padsumAmpMap = map;
padsumAmpMap(~pad_mask) = NaN;

% Use amplitude as weights
weights = padsumAmpMap;  %sum of amplitude weighted
%weights = map ./ (stdampMap.^2);   % charge weighted by inverse variance
weights(isnan(weights)) = 0;       % ignore NaNs

[xCentroid, yCentroid] = weightedAvg(Xgrid, Ygrid, weights);

figure;
subplot(1,2,1);

pad_xCut = xCenters <= pad_window_right & xCenters >= pad_window_left;
pad_projX = sum(padsumAmpMap, 1, 'omitnan');
bar(xCenters, pad_projX);
xlabel('X'); ylabel('Sum amplitude'); title('Projection on X (pad cut)');
hold on
fit_data = [];
fit_data(1,:) = xCenters(pad_xCut);
fit_data(2,:) = pad_projX(pad_xCut);
p0=[];
p0(1) = xCentroid;
p0(2) = 1;
p0(3) = max(pad_projX);
p0(4) = 0.01;
cmd='min; ret';
evalc("[p, err, chi] = fminuit('parabola4_minuit',p0,fit_data(:,1:end),'-b','-c',cmd);");
% store pad center x
pad_x_center = p(1);
pad_x_center_err = err(1);
% plot to see how parabolic fit looks like
plot(fit_data(1,:),parabola4_minuit(p, fit_data(1,:)),'LineWidth',2);

y_pos=get(gca,'ylim');
x_pos=get(gca,'xlim');
message = sprintf('Pad center x: %.2f \\pm %.2f\n', pad_x_center, pad_x_center_err);
message = [message sprintf('\\chi^2 / NDF: %.2f / %d', chi, length(fit_data(1,:)) - length(p))];
text(x_pos(1),0.95*y_pos(2),message);

grid on
hold off



subplot(1,2,2);

pad_yCut = yCenters <= pad_window_up & yCenters >= pad_window_bottom;
pad_projY = sum(padsumAmpMap, 2, 'omitnan');
bar(yCenters, pad_projY);
xlabel('Y'); ylabel('Sum amplitude'); title('Projection on Y (pad cut)');
hold on
fit_data = [];
fit_data(1,:) = yCenters(pad_yCut);
fit_data(2,:) = pad_projY(pad_yCut);
p0=[];
p0(1) = yCentroid;
p0(2) = 1;
p0(3) = max(pad_projY);
p0(4) = 0.01;
cmd='min; ret';
evalc("[p, err, chi] = fminuit('parabola4_minuit',p0,fit_data(:,1:end),'-b','-c',cmd);");

% store pad center y
pad_y_center = p(1);
pad_y_center_err = err(1);

% plot to see how parabolic fit looks like
plot(fit_data(1,:),parabola4_minuit(p, fit_data(1,:)),'LineWidth',2);

y_pos=get(gca,'ylim');
x_pos=get(gca,'xlim');
message = sprintf('Pad center y: %.2f \\pm %.2f\n', pad_y_center, pad_y_center_err);
message = [message sprintf('\\chi^2 / NDF: %.2f / %d', chi, length(fit_data(1,:)) - length(p))];
text(x_pos(1),0.95*y_pos(2),message);


grid on
hold off

pad_window_left = pad_x_center - tol*pad_window_length/2;
pad_window_right = pad_x_center + tol*pad_window_length/2;
pad_window_up = pad_y_center + tol*pad_window_length/2;
pad_window_bottom = pad_y_center - tol*pad_window_length/2;

pad_mask = Xgrid <= pad_window_right & Xgrid >= pad_window_left & Ygrid <= pad_window_up & Ygrid >= pad_window_bottom;

end

