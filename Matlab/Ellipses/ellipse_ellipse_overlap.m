function [overlap,area1,area2] = ellipse_ellipse_overlap(a1,b1,theta1,a2,b2,theta2,shouldShowFigure)
if nargin < 7
    shouldShowFigure = 0;
end

L = min(max(a1,b1),max(a2,b2));
n = 101;

% x=linspace(-L,L,n)';
% y=linspace(-L,L,n)';
% dx = x(2)-x(1);
% dy = y(2)-y(1);
% 
% [X,Y] = ndgrid(x,y);

x=linspace(-L,L,n);
y=linspace(-L,L,n).';
dx = x(2)-x(1);
dy = y(2)-y(1);
X = repmat(x,size(y));
Y = repmat(y,size(x));

XT = X*cos(theta1) + Y*sin(theta1);
YT = -X*sin(theta1) + Y*cos(theta1); 
ellipse1 = (XT/a1).^2 + (YT/b1).^2 < 1;
area1 = pi*a1*b1;

XT = X*cos(theta2) + Y*sin(theta2);
YT = -X*sin(theta2) + Y*cos(theta2); 
ellipse2 = (XT/a2).^2 + (YT/b2).^2 < 1;
area2 = pi*a2*b2;

overlap = sum(sum(ellipse1 & ellipse2))*dx*dy;

if shouldShowFigure == 1
    figure
    subplot(1,3,1)
    pcolor(X,Y,double(ellipse1)), shading flat
    subplot(1,3,2)
    pcolor(X,Y,double(ellipse2)), shading flat
    subplot(1,3,3)
    pcolor(X,Y,double(ellipse1 & ellipse2)), shading flat
end

end