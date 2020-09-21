function [overlap,area1,area2] = ellipse_ellipse_overlap_vector(a1,b1,theta1,a2,b2,theta2,shouldShowFigure)
if nargin < 7
    shouldShowFigure = 0;
end

n = 101;
alpha = linspace(0,2*pi-2*pi/n,n);
cosalpha = cos(alpha);
sinalpha = sin(alpha);
ellipse1 = polyshape(a1*cosalpha*cos(theta1) - b1*sinalpha*sin(theta1),a1*cosalpha*sin(theta1) + b1*sinalpha*cos(theta1));
ellipse2 = polyshape(a2*cosalpha*cos(theta2) - b2*sinalpha*sin(theta2),a2*cosalpha*sin(theta2) + b2*sinalpha*cos(theta2));

area1 = area(ellipse1);
area2 = area(ellipse2);
overlap = area(intersect(ellipse1,ellipse2));

if shouldShowFigure == 1
    figure
    plot(ellipse1), hold on
    plot(ellipse2)
    axis equal
end

end