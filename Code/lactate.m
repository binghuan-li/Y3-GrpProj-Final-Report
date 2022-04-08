Read= readmatrix('lactate.txt')
times = Read(:,1);
times = [times(36001:end)]
currents = Read(:,2);
currents = [currents(36001:end)]

%we want to cut all the extra stuff off the front, up to an hr which is
%3600*10= 36000th value

a= 3570;
b= 4880;
c= 6990;
d= 8360;
e= 9430;
f= 13600;
g= 14500;
h= 16950;
i= 19190;
j= 20230;
k= 24590;
l= 26190;
m= 27470;
n= 29120;
o= 30570;
p= 32000;
q= 33990;
r= 34840;

avg_current= zeros();

avg_current(1)= mean(currents(1:a));
avg_current(2)= mean(currents(a+1:b));
avg_current(3)= mean(currents(b+1:c));
avg_current(4)= mean(currents(c+1:d));
avg_current(5)= mean(currents(d+1:e));
avg_current(6)= mean(currents(e+1:f));
avg_current(7)= mean(currents(f+1:g));  
avg_current(8)= mean(currents(g+1:h));
avg_current(9)= mean(currents(h+1:i));
avg_current(10)= mean(currents(i+1:j));
avg_current(11)= mean(currents(j+1:k));
avg_current(12)= mean(currents(k+1:l));
avg_current(13)= mean(currents(l+1:m));
avg_current(14)= mean(currents(m+1:n));
avg_current(15)= mean(currents(n+1:o));
avg_current(16)= mean(currents(o+1:p));
avg_current(17)= mean(currents(p+1:q));
avg_current(18)= mean(currents(q+1:r));
avg_current(19)= mean(currents(r+1:end));


%%
%%clean single value step amperometry
trimmedtimes= [0,3570, 3571, 4880, 4881, 6990, 6991, 8360, 8361, 9430, 9431,13600, 13601, 14500, 14501, 16950, 16951, 19190, 19191, 20230, 20231,24590, 24591, 26190, 26191, 27470 , 27471, 29120, 29121, 30570, 30571, 32000, 32001,33990,33991, 34840, 34841,35493]

cur= [avg_current(1),avg_current(1), avg_current(2), avg_current(2), avg_current(3), avg_current(3), avg_current(4), avg_current(4), avg_current(5),avg_current(5), avg_current(6), avg_current(6), avg_current(7), avg_current(7), avg_current(8),avg_current(8), avg_current(9), avg_current(9), avg_current(10), avg_current(10), avg_current(11), avg_current(11), avg_current(12), avg_current(12), avg_current(13), avg_current(13), avg_current(14), avg_current(14), avg_current(15), avg_current(15), avg_current(16), avg_current(16), avg_current(17), avg_current(17), avg_current(18), avg_current(18), avg_current(19), avg_current(19)] 
 
 


plot(trimmedtimes,cur);
xlabel('time s','FontSize', 14)
ylabel('current A','FontSize', 14)
title('Lactate sensor amperometry', 'FontSize', 20)

%%calibration curve
conc=[0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.5,4.0,4.5]


plot(conc, avg_current,'b-d');
xlabel('concentration mM/L','FontSize', 14)
ylabel('Current A','FontSize', 14)
title('Lactate calibration curve','FontSize', 20)


%%line of best fit for linear range
x= [conc(1:10)]
y= [avg_current(1:10)]
coefficients = polyfit(x, y, 1);

xFit = linspace(min(x), max(x), 1000);

yFit = polyval(coefficients , xFit);

plot(x, y, 'b.', 'MarkerSize', 15); 
hold on; 
plot(xFit, yFit, 'r-', 'LineWidth', 2); 
xlabel('Lactate concentration mM/L','FontSize', 14);
ylabel('Current A','FontSize', 14);
title('Linear behaviour of lactate concentration,current behaviour','FontSize',20)

