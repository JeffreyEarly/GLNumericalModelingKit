function [overlap,area1,area2] = ellipse_ellipse_overlap_vector_custom(a1,b1,theta1,a2,b2,theta2,shouldShowFigure)
if nargin < 7
    shouldShowFigure = 0;
end

n = 101;
alpha = linspace(0,2*pi-2*pi/n,n).';
cosalpha = cos(alpha);
sinalpha = sin(alpha);
ellipse1.Vertices = [a1*cosalpha*cos(theta1) - b1*sinalpha*sin(theta1), a1*cosalpha*sin(theta1) + b1*sinalpha*cos(theta1)];
ellipse2.Vertices = [a2*cosalpha*cos(theta2) - b2*sinalpha*sin(theta2), a2*cosalpha*sin(theta2) + b2*sinalpha*cos(theta2)];

area1 = polyarea(ellipse1.Vertices(:,1),ellipse1.Vertices(:,2));
area2 = polyarea(ellipse2.Vertices(:,1),ellipse2.Vertices(:,2));

if area1 < area2
    % check if all of ellipse 1 is inside of ellipse 2
   isEllipse1PointsInsideEllipse2 = isInterior( ellipse2, ellipse1.Vertices(:,1),ellipse1.Vertices(:,2));
end

if area2 < area1
    % check if all of ellipse 1 is inside of ellipse 2
   isEllipse2PointsInsideEllipse1 = isInterior( ellipse1, ellipse2.Vertices(:,1),ellipse2.Vertices(:,2));
end 

if shouldShowFigure == 1
    figure
    plot(ellipse1.Vertices(:,1),ellipse1.Vertices(:,2)), hold on
    plot(ellipse2.Vertices(:,1),ellipse2.Vertices(:,2))
    axis equal
end

end

function isLeft = isLeft(polygon,i,x,y)
% using the notation above,
% -x3x4*y2y3 + x2x3*y3y4
isLeft = (polygon.dx(i)*(y-polygon.Vertices(i,2)) - (x-polygon.Vertices(i,1))*polygon.dy(i) );
end

function isinterior = isInterior(polygon, x, y)
% In my tests this was a 6-7x speedup over matlabs
% implementation.
windingNumber = zeros(size(x));
for i=1:(length(polygon.Vertices)-1)
    isless = polygon.Vertices(i,2) <= y;
    upwardCrossing = isless & polygon.Vertices(i+1,2) > y;
    if any(upwardCrossing)
        windingNumber(upwardCrossing) = windingNumber(upwardCrossing) + (IntegratorWithObstacles.isLeft(polygon,i,x(upwardCrossing),y(upwardCrossing)) > 0);
    end
    downwardCrossing = ~isless & polygon.Vertices(i+1,2) <= y;
    if any(downwardCrossing)
        windingNumber(downwardCrossing) = windingNumber(downwardCrossing) - (IntegratorWithObstacles.isLeft(polygon,i,x(downwardCrossing),y(downwardCrossing)) < 0);
    end
end
isinterior = abs(windingNumber) > 0;
end