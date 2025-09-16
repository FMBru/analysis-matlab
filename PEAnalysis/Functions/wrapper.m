% define array of runs you want to loop through, saving plots for
clear all;
close all;
runs = ['609';'610';'611';'612';'613';'614'];
i = 1;
    
for i=1:length(runs)
    %str_disp = sprintf('%s',runs(i));
    %disp(str_disp);
    current_run = runs(i,:);
    SPE_analysis;
end

