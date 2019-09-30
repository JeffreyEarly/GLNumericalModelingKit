theta1=25*pi/180;
a1=6; b1=3;

x=linspace(-10,10,101)';
y=linspace(-10,10,101)';
dx = x(2)-x(1);
dy = y(2)-y(1);

[X,Y] = ndgrid(x,y);
XT = X*cos(theta1) + Y*sin(theta1);
YT = -X*sin(theta1) + Y*cos(theta1); 
ellipse1 = sqrt((XT/a1).^2 + (YT/b1).^2) < 1;

area1 = pi*a1*b1;
area1_summed = sum(sum(ellipse1))*dx*dy;

theta2=120*pi/180;
a2=8; b2=2;

XT = X*cos(theta2) + Y*sin(theta2);
YT = -X*sin(theta2) + Y*cos(theta2); 
ellipse2 = sqrt((XT/a2).^2 + (YT/b2).^2) < 1;

figure
subplot(1,3,1)
pcolor(X,Y,double(ellipse1)), shading flat
subplot(1,3,2)
pcolor(X,Y,double(ellipse2)), shading flat
subplot(1,3,3)
pcolor(X,Y,double(ellipse1 & ellipse2)), shading flat

[overlap,area1,area2] = ellipse_ellipse_overlap(a1,b1,theta1,a2,b2,theta2);
