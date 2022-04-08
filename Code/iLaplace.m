function [g,yfit,config] = iLaplace(t,y,s,g0,alpha,mode)
%   HIGH LEVEL DESCRIPTION OF FUNCTION:
%   Inverse Laplace Transform with Regularized least squares 
%   algorithm to optimise

    g0 = g0(:); % [ ] initial guess of inverse laplace solution
    % Essentially the 'weight' of each exponential kernel component
    % at every instance of values in the s array
    s = s(:); % [s] range of lifetime values (time constant values)
    y = y(:); % [ ] chronoamperometric signal or current/nFA
    t = t(:); % [s] time
if strcmp(mode,'multiexp')
    plot_type = @semilogx;
    maxsearch = 20;
    options = optimset('MaxFunEvals',1e8,'MaxIter',1e8);
    constraints = {'g>0','zero_at_the_extremes'};
    R = zeros(length(g0)-2,length(g0));   % R is a matrix determining the form of the Regularizor.
    w = ones(size(y(:)));   % w is the array of fitting weigths for each value of y
end
    
[sM,tM] = meshgrid(s,t);
A = exp(-tM./sM); % Kernel in Fredholm integral equation
eps = 1e-09;      % threshold delta std to stop simulation
% Guess of g0
g0 = g0*sum(y)/sum(A*g0);  % Normalises g0 amplitude

% Plot Curve fit and g(s)
fh = gcf; clf(fh);
set(fh,'doublebuffer','on');
s1h = subplot(2,2,1); plot(t,y,'.',t,A*g0); title('Data and fitting curve'); axis tight
s2h = subplot(2,2,2); feval(plot_type,s,g0,'*-'); title('Initial distribution...'); axis tight
s3h = subplot(2,2,3); msdh = plot(0,0); title('Normalized msd');
s4h = subplot(2,2,4); plot(t,abs(y-A*g0)); title('Residuals'); axis tight

msd2plot = 0;
drawnow;
% Main cycle
ly = length(y);
oldssd = Inf;
tic

for k = 1:maxsearch,
    % msd: The mean square deviation. msd minimzed by fminsearch function
    % fminsearch uses unconstrained nonlinear minimization (Nelder-Mead)
    
    [g,msdfit] = fminsearch(@msd,g0,options,y,A,alpha,R,w,constraints);
    if ismember('zero_at_the_extremes',constraints)
        g(1) = 0;
        g(end) = 0;
    end    
    if ismember('g>0',constraints)
        g = abs(g);
    end

    g0 = g; % exponential kernel weights for the next step
    ssd = sqrt(msdfit/ly); % Sample Standard Deviation
    ssdStr = num2str(ssd);
    deltassd = oldssd-ssd; % Difference between "old ssd" and "current ssd"
    disp([int2str(k) ': ' ssdStr])
    oldssd = ssd;
    msd2plot(k) = msdfit/ly;
    plotdata(s1h,s2h,s3h,s4h,t,y,A,g,k,maxsearch,plot_type,s,msd2plot,msdh,ssdStr,deltassd);
    % Condition for the stabilization (end of cycles):
    % difference between "old ssd" and "current ssd" == 0
    if deltassd <= eps,
        disp(['Stabilization reached at step: ' int2str(k) '/' int2str(maxsearch)])
        break;
    end
end
disp(['Elapsed time: ' num2str(toc/60) ' min.'])
% Saving parameters and results
yfit = A*g; % fitting curve
config.t = t;
config.y = y;
config.yfit = yfit;
config.g0 = g0;
config.alpha = alpha;
config.R = R;
config.w = w;
config.maxsearch = maxsearch;
config.constraints = constraints;
config.date = datestr(now,30);
% Store g and s in a temporary file
% that can be loaded as a starting configuration
save iLaplace_temp g s config
disp('Temp output saved.')
set(gcf,'name','Fitting done')

% ### SUBS #####################################################
function out = msd(g,y,A,alpha,R,w,constraints)
% msd: The mean square deviation; this is the function
% that has to be minimized by fminsearch
% Constraints and any 'a priori' knowledge
if ismember('zero_at_the_extremes',constraints)
    g(1) = 0;
    g(end) = 0;
end
if ismember('g>0',constraints)
    g = abs(g); % must be g(i)>=0 for each i
end
r = diff(diff(g(1:end))); % second derivative of g
yfit = A*g;   %% Series of exponentials
% Sum of weighted square residuals
VAR = sum(w.*(y-yfit).^2);
% Regularizor
REG = alpha^2 * sum((r-R*g).^2);
% Output to be minimized
out = VAR+REG;

% ### SUBS #####################################################
function plotdata(s1h,s2h,s3h,s4h,t,y,A,g,k,maxsearch,plot_type,s,msd2plot,msdh,ssdStr,deltassd)
% For the "real-time" plots
axes(s1h);
plot(t,y,'.',t,A*g); title('Data')
xlabel('Time (s)');
axis tight
axes(s2h);
feval(plot_type,s,g,'o-');
title('Relaxation times distribution g(s)');
xlabel('s');
axis tight
axes(s3h);
title(['ssd: ' ssdStr '; \Deltassd: ' num2str(deltassd)])
ylabel('Sample Standard Deviation')
xlabel('Step')
set(msdh,'xdata',1:k,'ydata',msd2plot);
axis tight
axes(s4h);
plot(t,abs(y-A*g)/length(y),'o-');
title('Normalized residuals')
xlabel('Time (s)');
axis tight
set(gcf,'name',['Step:' int2str(k) '/' int2str(maxsearch)])
drawnow