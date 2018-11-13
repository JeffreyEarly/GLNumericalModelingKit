function [p, p_err] = LinearLeastSquaresFit(X,Y,shouldSetInterceptToZero)
% Returns the least squares fit with standard error.
%
% If shouldSetInterceptToZero is empty or zero, the data will be fit to the
% function p(2)*X+p(1). The standard error for the coefficients are
% returned in the same order.
%
% If shouldSetInterceptToZero=1, the data will be fit to p(1)*X.
if nargin < 3
    shouldSetInterceptToZero = 0;
end
if shouldSetInterceptToZero == 0
    %% Calculation of Standard Error With Intercept
    n = length(X);                                      %  Number of observations
    XBar=mean(X);                                       %  Calculates mean of X
    YBar=mean(Y);                                       %  Calculates mean of Y
    Sxx=sum((X-XBar).^2);
    Sxy=sum((X-XBar).*(Y-YBar));
    Slope = Sxy/Sxx;                                    %  Calculates Slope
    Intercept= YBar-Slope*XBar;                         %  Calculates Intercept
    yfit=Intercept + X*Slope;                           %  Fitted response values based on the slope
    r = Y - yfit;                                       %  r is the residuals, which is the observed minus fitted values
    SSE = sum(r.^2);                                    %  SSE is the sum of squared errors
    MSE=SSE/(n-2);                                      %  Mean Squared Error
    p_err =[                                 %  Standard Error of the regression coefficients
        sqrt(MSE*sum(X.^2)/(n*Sxx));                    %  Standard Error of the intercept coefficient
        sqrt(MSE/Sxx)];                                  %  Standard Error of the slope coefficient
    p = [Intercept; Slope];
elseif shouldSetInterceptToZero == 1
    %% Calculation of Standard Error Without Intercept
    % https://www.mathworks.com/matlabcentral/answers/373940-how-does-matlab-calculate-standard-error-in-fitlm
    n = length(X);                                       %  Number of observations
    Slope = sum(X.*Y)/sum(X.^2);                         %  Calculates slope
    yfit=X*Slope;                                        %  Fitted response values based on the slope
    r = Y - yfit;                                        %  r is the residuals, which is the observed minus fitted values
    SSE = sum(r.^2);                                     %  SSE is the sum of squared errors
    MSE=SSE/(n-1);                                       %  Mean Squared Error
    SE=sqrt(MSE/sum(X.^2));
    
    p = Slope;
    p_err = SE;
else
    error('invalid choice') 
end
    
