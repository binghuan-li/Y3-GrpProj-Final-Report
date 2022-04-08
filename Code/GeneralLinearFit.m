function [FitOut]=GeneralLinearFit(x,y)

%% fits data to y=b0+b1*x and returns 95% ME and confidence bands

[xData, yData] = prepareCurveData( x(:), y(:) );
ft = fittype( 'poly1' );
opts = fitoptions( ft );
[fitresult] = fit( xData, yData, ft, opts );   % Fit model to data. fitresult shows the output of the fit

[p]=coeffvalues(fitresult);         % get coefficient values
b1=p(1);                            % extract the slope 
b0=p(2);                            % extract the intercept

CI=confint(fitresult);              % first row is lower bound of CI of b1 then b0, second row is upper bound of b1 then b0
ME1=CI(2,1)-b1;                     % Calculate Margin of Error (half width of CI)
ME0=CI(2,2)-b0;                     % Calculate Margin of Error (half width of CI)

xfit=linspace(min(x),max(x),1000);                  % create an array of 1000 values beteween the minimum and maximum value of x
yfit=feval(fitresult,xfit);                         % calculate the model for the new x values
cbands=predint(fitresult,xfit,0.95,'functional');   % calculate the 95% confidence bands for the fit
pbands=predint(fitresult,xfit,0.95,'observation');  % calculate the 95% confidence bands for the fit

% Built the structure with the output parameters of interest
FitOut.b0=b0;
FitOut.b1=b1;
FitOut.ME0=ME0;
FitOut.ME1=ME1;
FitOut.xfit=xfit;
FitOut.yfit=yfit;
FitOut.cbands=cbands;
FitOut.pbands=pbands;
