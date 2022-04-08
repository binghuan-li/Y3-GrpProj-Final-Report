clear all
close all
clc

%% Code to Extract Electron transfer rate constants from Chronoamperometric Data
% Example Equation: I/nFA = (A_B.*exp(-k_AT.*(t))+A_U.*exp(-k_A.*(t))) 
k_A_v = [200 250 150 100];
k_AT_v = [400 350 450 500];
for i = 1:size(k_A_v,2)
    % Parameters
    A_O = 1;            % [ ] Total Aptamer Concentration
    A_B = 0.5*A_O;      % [ ] Bound Aptamer Concentration
    A_U = 0.5*A_O;      % [ ] Unbound Aptamer Concentration
    k_A = k_A_v(i);            % [s^-1] Electron rate transfer constant unbound state
    k_AT = k_AT_v(i);           % [s^-1] Electron rate transfer constant bound state
    sol_start = log10(1/800); % 
    sol_end = log10(1/80);
    %sol_start = -3;
    %sol_end = -1;
    % s space (range of possible lifetime values):
    s0 = logspace(sol_start,sol_end,150)';

    % Define time array
    h = 0.00001; % step size ASK --- WHAT IS THE TIME STEP OF MEASUREMENTS?
    t = (0:h:0.1)'; % time array
    snr = 28; % Signal to Noise Ratio?
    [sM,tM] = meshgrid(s0,t); % ignore 
    % Data:
    y = (A_B.*exp(-k_AT.*(t))+A_U.*exp(-k_A.*(t)));
    y_n = awgn(y,snr,'measured'); % the signal level of y is computed to determine appropriate noise level based on value of snr
    % Regularisation Parameter
    alpha = 1;

    % Guess s space and g function
    s = logspace(sol_start,sol_end,150)';
    g = ones(size(s));
    [g,yfit,cfg] = iLaplace(t,y_n,s,g,alpha,'multiexp');

    figure (2)
    subplot(1,2,1)
    plot(t,y_n,'b.','MarkerSize',10);
    hold on
    plot(t,yfit,'-r','MarkerSize',4.0,'LineWidth',1.5);
    hold on
    title('Raw data and Fit','interpreter','latex')
    ylabel('Concentration of Biomarker []','interpreter','latex', 'Fontsize',16,'interpreter','latex')
    xlabel('Time [s]','Fontsize', 16,'interpreter','latex')
    set(gca,'color','w')
    legend('Curve fit','Raw Data','interpreter','latex','Fontsize',12)
    % f1 = fit(1./s(1:35),g(1:35)./max(g),'gauss2')
    % f2 = fit(1./s(60:95),g(60:95)./max(g),'gauss2')
    subplot(1,2,2)
    plot(1./s,g/max(g),'-*'); title('Actual Decay rates, Algorithm Decay Rates','interpreter','latex','FontSize',15);
    hold on
    plot(1./s,g/max(g),'-*'); title('Actual Decay rates, Algorithm Decay Rates','interpreter','latex','FontSize',15);


    hold on

    xline(k_A,'--r', 'LineWidth', 1)
    xline(k_AT, '--r','LineWidth', 1)
    legend('Calculated rate constants','$\tau_{1} = 200s^{-1}$','$\tau_{2} = 400s^{-1}$','interpreter','latex','Fontsize',12)
    ylabel('Normalised Intensity','Fontsize',16,'interpreter','latex')
    xlabel('$Rate Constant [s^{-1}]$','Fontsize',16,'interpreter','latex')
    set(gca,'color','w')
    
    fname = sprintf('Alpha_testing%d_%d150points.mat',k_A_v(i),k_AT_v(i));
    save(fname)
end
