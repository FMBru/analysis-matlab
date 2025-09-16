% Author: M.L. Velazquez Fernandez
% 
% Script created to plot the SPE mean vs the cathode V
% and then fit an exponential to it.

clear;
clc;

% Define the data points
x = [275, 285, 295, 305, 315] % Cathode voltage
y = [0.01364, 0.01825, 0.02396, 0.03150, 0.04113]; % SPE mean
y_errors = [0.00042, 0.00006, 0.00009, 0.00014, 0.00023]

% Plot the datapoints
figure;
errorbar(rdivide(x,500), y, y_errors, 'vertical')
scatter(rdivide(x,500), y,'filled', 'o', 'black');
set(gca, 'yscale', 'log')
hold on;

% Fit and plot the exponential to the data
f = fit(rdivide(x,500)', y', 'exp1');
coeffs = coeffvalues(f);
a = coeffs(1);
b = coeffs(2);
eq = sprintf('y = %.4e * exp(%.4e * x)', a, b)
plot(f, '-- black')
set(gca, 'yscale', 'log')
% Elec. mesh with anode voltage fix at 300V
x_mesh = [510, 520, 530, 540]
y_mesh = [0.01634, 0.02127, 0.02655, 0.03365]

scatter((275 ./ x_mesh), y_mesh, 'filled', 'o', "MarkerFaceColor","r");

f1 = fit((275 ./ x_mesh)', y_mesh', 'exp1');
coeffs1 = coeffvalues(f1);
a1 = coeffs1(1);
b1 = coeffs1(2);
eq1 = sprintf('y = %.4e * exp(%.4e * x)', a1, b1);

plot(f1, '-- red');
set(gca, 'yscale', 'log')
% Elec. mesh with anode voltage fix at 300V
x_mesh1 = [500, 510, 520]
y_mesh1 = [0.02601, 0.03256, 0.03948]

scatter(rdivide(300, x_mesh1), y_mesh1, 'filled', 'o', "MarkerFaceColor","b");

f2 = fit(rdivide(300, x_mesh1)', y_mesh1', 'exp1');
coeffs2 = coeffvalues(f2);
a2 = coeffs2(1);
b2 = coeffs2(2);
eq2 = sprintf('y = %.4e * exp(%.4e * x)', a2, b2);

plot(f2, '-- blue');



%
xlabel('Voltage ratio [Anode/Cathode]');
ylabel('SPE mean');

legend('Elec. Mesh (500C)', ['Exponential fit: ', eq], 'Elec. Mesh (275A)', ['Exponential fit: ', eq2], 'Elec. Mesh (300A)', ['Exponential fit: ', eq2]); 
hold off;



