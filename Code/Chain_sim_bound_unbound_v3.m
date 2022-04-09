%% Simulating Aptamers using Freely-jointed chain model
clear all
close all
clc

% Persistence length of single DNA = 1-4nm
% Bernard, Tinland (1997). "Persistence Length of Single-Stranded DNA". Macromolecules. 30 (19): 5763
% Kuhn length is considered to be twice the persistence length
% 

% probability distribution of distance of redox reporter from the 
N = 5E6; % [ ] number of simulations
phi_max = 180;
theta_max = 360;
l = 0.34e-9; % [m] = 0.34nm length of one base
pers_length = 2e-09;
r = 2*pers_length; % [nm] kuhn length
bases_bound = 50;
bases_unbound = 100;
DNA_length_bound  = ceil((bases_bound*l)/r); % [ ] number of kuhn lengths for unbound Aptamer
DNA_length= ceil((bases_unbound*l)/r);% [ ] number of kuhn lengths for unbound Aptamer
theta_v = zeros(1,DNA_length);
phi_v = zeros(1,DNA_length);
length_check = zeros(1,DNA_length);
x_v = zeros(1,DNA_length+1);
y_v = zeros(1,DNA_length+1);
z_v = zeros(1,DNA_length+1);
x_P_b = zeros(1,DNA_length_bound+1);
y_P_b = zeros(1,DNA_length_bound+1);
z_P_b = zeros(1,DNA_length_bound+1);
origin = [0 0 0];
plot_cont = 1; 
x_temp = 0;
y_temp = 0;
z_temp = 0;
apt_count = 0;
theta_v(1) = 0;
phi_v(1) = 0;
x_v(1) = 0;
y_v(1) = 0;
z_v(1) = 0;
x_P_b(1) = 0;
y_P_b(1) = 0;
z_P_b(1) = 0;
theta_v2(1) = 0;
phi_v2(1) = 0;
counter = 0;
sim_no = 1;

while sim_no <= N
    for k = 2:(DNA_length+1)
        theta_rand = rand*(theta_max);
        phi_rand = rand*(phi_max);
        x_temp = x_v(k-1) + r*sind(phi_rand)*cosd(theta_rand);
        y_temp = y_v(k-1) + r*sind(phi_rand)*sind(theta_rand);
        z_temp = z_v(k-1) + r*cosd(phi_rand);
        z_sum(k-1) = z_temp;
        if z_temp < 0
            z_temp_2 = z_temp;
            z_temp = abs(z_temp_2);
            phi_rand = acosd((z_temp-z_v(k-1))/r);
            x_temp = x_v(k-1) + r*sind(phi_rand)*cosd(theta_rand);
            y_temp = y_v(k-1) + r*sind(phi_rand)*sind(theta_rand);
        end
        theta_v(k-1) = theta_rand; % Storing the choice of angle
        phi_v(k-1) = phi_rand;
        x_v(k) = x_temp;
        y_v(k) = y_temp;
        z_v(k) = z_temp;
        length_check(k-1) = ((x_temp-x_v(k-1))^2+(y_temp-y_v(k-1))^2+(z_temp-z_v(k-1))^2)^(1/2);
    end
        apt_count = apt_count+1;
        z_distance(apt_count) = z_temp;
%         figure (1)
%         plot3(x_v,y_v,z_v,'MarkerSize',2,'MarkerFaceColor','k','LineWidth',1.5)
%         hold on
%         plot3(x_v(1,DNA_length+1),y_v(1,DNA_length+1),z_v(1,DNA_length+1),'bo','MarkerSize',10,'MarkerFaceColor','b')
%         hold on
        counter = counter+1;
        x_temp = 0;
        y_temp = 0;
        z_temp = 0;
%         x_mat(:,sim_no) = x_v;
%         y_mat(:,sim_no) = y_v;
%         z_mat(:,sim_no) = z_v;
        length_check_mat(:,sim_no) = length_check;
        x_v = zeros(1,DNA_length+1);
        y_v = zeros(1,DNA_length+1);
        z_v = zeros(1,DNA_length+1);
        plot_cont = 1;
        sim_no = sim_no+1;
end

z_avg = mean(z_distance)
B = sort(z_distance);
figure (2)
plot(z_distance);


figure (3)
numIntervals = 100;
intervalWidth = (max(B) - min(B))/numIntervals;
x = 0:intervalWidth:1.5e-08;
ncount = histc(B,x);
relativefreq = ncount/length(B);
bar(x-intervalWidth/2, relativefreq,1)
xlim([min(x) max(x)])
% set(gca, 'xtick', x)
hold on

%%
    apt_count = 0;
    x_temp = 0;
    y_temp = 0;
    z_temp = 0;
    N = 5E6; % [ ] number of simulations
    plot_count = 1;
    counter_b = 0;
    sim_no = 1;
    figure (4)

while sim_no <= N
    for k = 2:(DNA_length_bound+1)
        theta_rand = rand*(theta_max);
        phi_rand = rand*(phi_max);
        x_temp = x_P_b(k-1) + r*sind(phi_rand)*cosd(theta_rand);
        y_temp = y_P_b(k-1) + r*sind(phi_rand)*sind(theta_rand);
        z_temp = z_P_b(k-1) + r*cosd(phi_rand);
        z_sum(k-1) = z_temp;
        if z_temp < 0
            z_temp_2 = z_temp;
            z_temp = abs(z_temp_2);
            phi_rand = acosd((z_temp-z_v(k-1))/r);
            x_temp = x_P_b(k-1) + r*sind(phi_rand)*cosd(theta_rand);
            y_temp = y_P_b(k-1) + r*sind(phi_rand)*sind(theta_rand);
        end
        theta_v2(k-1) = theta_rand; % Storing the choice of angle
        phi_v2(k-1) = phi_rand;
        x_P_b(k) = x_temp;
        y_P_b(k) = y_temp;
        z_P_b(k) = z_temp;
        length_check(k-1) = ((x_temp-x_P_b(k-1))^2+(y_temp-y_P_b(k-1))^2+(z_temp-z_P_b(k-1))^2)^(1/2);
    end
        apt_count = apt_count+1;
        z_distance_bound(apt_count) = z_temp;
%         figure (4)
%         plot3(x_P_b,y_P_b,z_P_b,'MarkerSize',2,'MarkerFaceColor','k','LineWidth',1.5)
%         hold on
%         plot3(x_P_b(1,DNA_length_bound+1),y_P_b(1,DNA_length_bound+1),z_P_b(1,DNA_length_bound+1),'bo','MarkerSize',10,'MarkerFaceColor','b')
%         hold on
        counter = counter+1;
        x_temp = 0;
        y_temp = 0;
        z_temp = 0;
        length_check_mat(:,sim_no) = length_check;
%         x_mat(:,sim_no) = x_v;
%         y_mat(:,sim_no) = y_v;
%         z_mat(:,sim_no) = z_v;
        x_P_b = zeros(1,DNA_length_bound+1);
        y_P_b = zeros(1,DNA_length_bound+1);
        z_P_b = zeros(1,DNA_length_bound+1);
        plot_cont = 1;
        sim_no = sim_no+1;
end

    z_avg_2 = mean(z_distance_bound)
    B_bound = sort(z_distance_bound);
    figure (2)
    plot(z_distance_bound)


figure (3)
numIntervals = 100;
% intervalWidth = (max(B_bound) - min(B_bound))/numIntervals;
x = 0:intervalWidth:1.5e-08;
ncount = histc(B_bound,x);
relativefreq_bound = ncount/length(B_bound);
peter = bar(x-intervalWidth/2, relativefreq_bound,1,'FaceAlpha',0.5)
xlim([min(x) max(x)])
xlabel('length [m]', 'FontSize', 5)
% set(gca, 'xtick', x)
ax = gca;
ax.FontSize = 10; 
legend('30 bases Freely Jointed Chains (FJCs)', '10 bases FJCs')
% bar(z_distance')
x_m = x*(10^9);

prob_unbound = sum(relativefreq(1:12));
prob_bound = sum(relativefreq_bound(1:12));



