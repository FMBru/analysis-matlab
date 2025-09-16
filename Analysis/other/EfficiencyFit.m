close all;
clear all;

x = [1.9, 3.5, 3.4, 2.0, 5.25, 20.2];
y = [81.1, 91.8, 95.6, 83.8, 99.4, 100.0];

%x = [1, 4, 9, 16, 25]
%y = [1, 2, 3, 4, 5]

sqrtFit = fittype( @(a,b,x) a*x.^(0.5)+b);

[f,gof] = fit(x(:),y(:),sqrtFit)

P = [f.a, f.b];

figure
hold on
plot(f,x,y,"o");
grid on;
xlabel('NPE')
ylabel('Efficiency [%]')
xlim([0 25]);
ylim([80 100]);

lgd = legend;
lgd.Location = "southeast";

caption = sprintf([' f(x) = %.2f * x^{0.5} + %.2f  \n R^2 = %.2f '], P(1), P(2), gof.rsquare);
text(4, 92, caption);

