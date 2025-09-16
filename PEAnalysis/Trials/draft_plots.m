% Author: M.L. Velazquez Fernandez
% 
% Script created to plot the SPE mean vs the cathode V
% and then fit an exponential to it.

clear;
clc;

% Define the data points (thin mesh)
x = [530,550,570,590,600,610,620] % Cathode voltage
y = [0.00394, 0.00597, 0.00962, 0.01587, 0.02096, 0.02711, 0.03675]; % SPE mean
y_errors = [0.00042, 0.00006, 0.00009, 0.00014, 0.00023]

% Plot the datapoints
figure;
%errorbar(rdivide(x,500), y, y_errors, 'vertical')
scatter(x, y,'filled', 'o', 'black');
%set(gca, 'yscale', 'log')
hold on;

% Fit and plot the exponential to the data
f = fit(x', y', 'exp1');
coeffs = coeffvalues(f);
a = coeffs(1);
b = coeffs(2);
eq = sprintf('y = %.4e * exp(%.4e * x)', a, b)
plot(f, '-- black')
%set(gca, 'yscale', 'log')

% Define the data points (std mesh)
x1 = [530,550,570,590,610,620] % Cathode voltage
y1 = [0.00275, 0.0043, 0.00738, 0.01251, 0.02252, 0.02824]; % SPE mean

x2 = [540,560,580,600] % Cathode voltage
y2 = [0.00297, 0.00609, 0.01006, 0.01666]; % SPE mean

% Plot the datapoints
%errorbar(rdivide(x,500), y, y_errors, 'vertical')
scatter(x1, y1,'filled', 'o', 'red'); % 20/06/24
scatter(x2, y2,'filled', '^', 'red'); % 14/06/24
%set(gca, 'yscale', 'log')

% Fit and plot the exponential to the data
f1 = fit(x1', y1', 'exp1');
coeffs1 = coeffvalues(f1);
a1 = coeffs1(1);
b1 = coeffs1(2);
eq1 = sprintf('y = %.4e * exp(%.4e * x)', a1, b1)
plot(f1, '- red')
%set(gca, 'yscale', 'log')




% % Define the data points
% x = [10, 30, 70, 100, 150 , 200] % Cathode voltage
% y = [0.00578, 0.00585, 0.00605, 0.00613, 0.00683, 0.01001]; % SPE mean
% 
% % Plot the datapoints
% figure;
% %errorbar(rdivide(x,500), y, y_errors, 'vertical')
% scatter(x, y,'filled', 'o', 'black');
% %set(gca, 'yscale', 'log')
% hold on;
% 
% % Fit and plot the exponential to the data
% f = fit(x', y', 'exp2');
% coeffs = coeffvalues(f);
% a = coeffs(1);
% b = coeffs(2);
% eq = sprintf('y = %.4e * exp(%.4e * x)', a, b)
% plot(f, '-- black')
% %set(gca, 'yscale', 'log')
% 
% 
% % Define the data points
% x1 = [10, 30, 50, 70, 100, 150 , 200] % Cathode voltage
% y1 = [0.00941, 0.00958, 0.00962, 0.00967, 0.00982, 0.01149, 0.01715]; % SPE mean
% 
% scatter(x1, y1,'filled', 'o', 'red');
% %set(gca, 'yscale', 'log')
% hold on;
% 
% % Fit and plot the exponential to the data
% f1 = fit(x1', y1', 'exp2');
% coeffs1 = coeffvalues(f1);
% a1 = coeffs1(1);
% b1 = coeffs1(2);
% eq1 = sprintf('y = %.4e * exp(%.4e * x)', a1, b1)
% plot(f1, '-- red')
% %set(gca, 'yscale', 'log')
% 
% 
% % Define the data points (STD mesh 550A)
% x2 = [ 30, 70, 100, 150 , 200] % Cathode voltage
% y2 = [0.00519, 0.00534, 0.00573, 0.00773, 0.01591]; % SPE mean
% 
% scatter(x2, y2,'filled', '^', 'black');
% %set(gca, 'yscale', 'log')
% hold on;
% 
% % Fit and plot the exponential to the data
% f2 = fit(x2', y2', 'exp2');
% coeffs2 = coeffvalues(f2);
% a2 = coeffs2(1);
% b2 = coeffs2(2);
% eq2 = sprintf('y = %.4e * exp(%.4e * x)', a2, b2)
% plot(f2, '- black')
% %set(gca, 'yscale', 'log')
% 
% % Define the data points (STD mesh 570A)
% x3 = [10, 30, 70, 100, 150 , 200] % Cathode voltage
% y3 = [0.00725, 0.00762, 0.00816, 0.00885, 0.01259, 0.02914]; % SPE mean
% 
% scatter(x3, y3,'filled', '^', 'red');
% %set(gca, 'yscale', 'log')
% hold on;
% 
% % Fit and plot the exponential to the data
% f3 = fit(x3', y3', 'exp2');
% coeffs3 = coeffvalues(f3);
% a3 = coeffs3(1);
% b3 = coeffs3(2);
% eq3 = sprintf('y = %.4e * exp(%.4e * x)', a3, b3)
% plot(f3, '- red')
% %set(gca, 'yscale', 'log')
% 
xlabel('Cathode Voltage [V]');
ylabel('SPE mean');

%legend('Thin. Mesh (550V A)', ['Exponential fit: ', eq], 'Thin. Mesh (570V A)', ['Exponential fit: ', eq1], 'STD Mesh (550V A)', ['Exponential fit: ', eq2], 'STD Mesh (570V A)', ['Exponential fit: ', eq3]); 

legend('Thin. Mesh (50V C)', ['Exponential fit: ', eq], 'STD Mesh (50V C;  20/06/24)', 'STD Mesh (50V C;  14/06/24)', ['Exponential fit: ', eq1])
hold off;

