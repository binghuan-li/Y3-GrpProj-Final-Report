%% coulometry calculation 
% Binghuan Li - Group 5, Biomedical Engineering Year-3 Project 
% Dept. of Bioengineering, Imperial College London
% ============RELEASE NOTES==========================
% INITIAL Ver0: 22-03-2022
% UPDATE NOTES: 
% Ver1 - 03-04-2022 fixed linear fitting range; 
%                   fixed legend mismatch
% Ver2 - 04-04-2022 implemented concentration calc; 
%                   refined linear fitting range;
% ====================================================
clc;
close all;
clear all;

%% Construct Anson's Equation, all in SI units
F = 96485.3365; % Faraday constant 
D = 1.75E-6;    % diffusion coefficient of H2O2, https://doi.org/10.1016/j.sbsr.2018.10.001
n = 2;          % electron transfer No.
r = 1*(10^-3);  % radius of eletrode - 1mm
A = r.*r.*pi;   % area of electrode

counter = 1;
figure(1);
for c = [0, 50, 100, 150, 200, 250, 300] % concentration range
    % read csv files
    file_name = append(int2str(c),'.xlsx');
    amp_data = xlsread(file_name);
    
    % get corrected time step due to the edge effect
    t = amp_data(:,1);
    sliced_t = t(1:6000,1); % only need first 6000 points
    current = amp_data(:,2);
    current_sliced = (current(1:6000,1).*(1E-3)); % only need first 6000 points
    edge_corrected_t = sqrt(sliced_t)+(((1.92.*sqrt(D))./r).*sliced_t);
    
    % integration
    Q = cumtrapz(sliced_t,current_sliced);
    
    % plot figures 
    plt(counter) = plot(edge_corrected_t,Q, 'linewidth', 2);
    hold on;
    
    % fitting
    Q_linear = Q(1000:1500,1); % linear fitting based on the grad between pt 1000-1500
    t_linear =  edge_corrected_t(1000:1500,1);
    ft = fittype('a*x+b');
    f{counter} = fit(t_linear,Q_linear,ft);
    plot(f{counter});
    disp(f{counter});
    
    % calculate concentration in uM
    validated_conc{counter} = (((f{counter}(2)-f{counter}(1))./(2-1))./((n.*F.*A*sqrt(D))./(sqrt(pi))))*(1E-3);
    
    counter = counter + 1;
end
hold off;
xlabel("Edge corrected time");
ylabel("Q");
title("Hydrogen Peroxide Coulometry Plot")
legend([plt(1) plt(2) plt(3) plt(4) plt(5) plt(6) plt(7)], {"0uM", "50uM", "100uM", "150uM", "200uM", "250uM", "300uM"});

disp(validated_conc);

