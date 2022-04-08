%% Electron transfer occurs when the aptamer is within a distance threshold
clear all
close all
clc

load('Aptamer_Gaussian_Lengths_Data.mat')
% Gaussian Curve Fits
g_50=A1.*exp((-(y-mu1).^2)./(2*sigma1^2));
g_100=A2.*exp((-(y-mu2).^2)./(2*sigma2^2));
g_200=A3.*exp((-(y-mu3).^2)./(2*sigma3^2));
% Plots
figure (1)
plot(y,g_100);
hold on
plot(y,g_50);
% Assume DNA chains are non-conducting. Electron transfer happens or
% doesn't. Use a threshold distance to quantify number of electrons within
% that threshold
% Parameters
dist_thresh = 0.5*10^-9; % [m]
% x = 0:0.1:20;
nApt(1) = 10*10^-9; %total number of aptamers 
bRatio = 1;
n=2; %electrons per redox probe 
F= 96485; %faradays constant 
A = pi*1.5*10^-3; %area of electrode 
T = 30;
k_ET = 99999; % infinite rate constant for aptamers within distance threshold
Charge = 1.602176634E-19; % Charge in a single electron
deltat = 0.01;
time = zeros(T/deltat,1);
i = zeros(T/deltat,1);

% Probabilities
p_unbound = zeros(length(x),1);
p_bound = zeros(length(x),1);
p_bound_tot = 0;
p_unbound_tot = 0;


for i = 1:size(x,2)
    if x(i) <= 1.5
        p_bound_thresh_50(i) = relativefreq_50(i);
        p_unbound_thresh_100(i) = relativefreq_100(i);
        p_unbound_thresh_200(i) = relativefreq_200(i);

    else
        break
    end
end

load('Aptamer_Gaussian_Lengths_Data_no2.mat')
for i = 1:size(x,2)
    if x(i) <= 1.5
        p_bound_thresh_10(i) = relativefreq_10(i);
        p_bound_thresh_20(i) = relativefreq_20(i);
        p_bound_thresh_30(i) = relativefreq_30(i);
        p_bound_thresh_40(i) = relativefreq_40(i);
    else
        break
    end
end

load('Aptamer_Gaussian_Lengths_Data_no3.mat')
for i = 1:size(x,2)
    if x(i) <= 1.5
        p_bound_thresh_80(i) = relativefreq_80(i);
        p_bound_thresh_150(i) = relativefreq_150(i);
    else
        break
    end
end

p_bound_tot_10 = sum(p_bound_thresh_10);

p_bound_tot_20 = sum(p_bound_thresh_20);

p_bound_tot_30 = sum(p_bound_thresh_30);

p_bound_tot_40 = sum(p_bound_thresh_40);

p_bound_tot_50 = sum(p_bound_thresh_50);

p_bound_tot_80 = sum(p_bound_thresh_80);

p_unbound_tot_150 = sum(p_bound_thresh_150);

p_unbound_tot_100 = sum(p_unbound_thresh_100);
p_unbound_tot_200 = sum(p_unbound_thresh_200);

lengths = [10 20 30 40 50 80 100 150 200];
probabilities = [p_bound_tot_10 p_bound_tot_20 p_bound_tot_30 p_bound_tot_40 p_bound_tot_50 p_bound_tot_80 p_unbound_tot_100 p_unbound_tot_150 p_unbound_tot_200]

figure (500)
plot(lengths, probabilities,'-s', 'LineWidth',2)
hold on
set(gcf,'color','w')
xlabel('Number of Bases')
ylabel('Probability of Collision [ ]')
title('Dependance of Probability of Collision (L < 1.5nm) on Chain Length')
grid on
box on


figure (501)
loglog(lengths, probabilities,'-s', 'LineWidth',2)
hold on
set(gcf,'color','w')
xlabel('Number of Bases')
ylabel('Probability of Collision [ ]')
title('Dependance of Probability of Collision (L < 1.5nm) on Chain Length log-log plot')
hold on

f = GeneralLinearFit(log10(lengths(4:9))',log10(probabilities(4:9))');
    p1 = f.b0;
    p2 = f.b1;
    c_bands_N = f.cbands;
    xfit_3 = f.xfit;
    yfit_3 = f.yfit;
    figure (501)
    hold on
    loglog(10.^(xfit_3),10.^(yfit_3),'-g','LineWidth', 1.4)
    hold on
    loglog(10.^(xfit_3)',10.^(c_bands_N(:,2)),'--r','LineWidth', 1.2);
    hold on
    loglog(10.^(xfit_3)',10.^(c_bands_N(:,1)), '--r','LineWidth', 1.2);
    set(gca,'fontsize',18,'TickLabelInterpreter','latex')
%     legend('Raw Data','Excluded in Curve fit',sprintf('Curve fit: log10(P(L<1.5)) = %f $\times$ N + %f', p2, p1),'FontSize',15, 'Interpreter','latex');
    legend('Raw Data','Curve fit: $log10(P(L<1.5)) = -0.4272 \times N - 0.0056$','FontSize',15, 'Interpreter','latex');

disp(p1)
disp(p2)


% Assuming Electron transfer only occurs below a threshold distance (1nm)
% Electron transfer is Nernstian within this Threshold (effectively
% instantaneous). k[C] . The rates equal to concentration [AT]. effective rate
% constant. Count number of Ends of Aptamers within the threshold (1nm).
% Rate constant there of electron transfer is infinite, therefore Nernst
% equation applies. Effective rate constant is therefore proportional to total number of
% particles within that threshold distance. The probability of distances
% will therefore come in as the concentration term in the multiexponential equation
% Current is the rate of transfer of electrons

