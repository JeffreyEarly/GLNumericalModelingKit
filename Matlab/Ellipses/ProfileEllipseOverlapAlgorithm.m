theta1=0*pi/180;
a1=6; b1=3;

theta2=45*pi/180;
a2=8; b2=2;

for i=1:1000
    [overlap,area1,area2] = ellipse_ellipse_overlap_vector(a1,b1,theta1,a2,b2,theta2,0);
end
