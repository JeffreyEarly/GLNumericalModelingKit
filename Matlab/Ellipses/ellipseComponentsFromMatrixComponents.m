function [semiMajor,semiMinor,theta] = ellipseComponentsFromMatrixComponents(Mxx,Myy,Mxy)
%ellipseComponentsFromMatrixComponents Give matrix components Mxx,Myy,Mxy,
% this returns the equivalent representation as a rotated ellipse in the
% form of semi-major axis, semi-minor axis, and rotation angle.
descriminant = sqrt((Mxx-Myy).*(Mxx-Myy) + (4*Mxy.*Mxy));
D2_plus = (Mxx+Myy+descriminant)/2;
D2_minus = (Mxx+Myy-descriminant)/2;
semiMajor = sqrt(D2_plus);
semiMinor = sqrt(D2_minus);
theta = atan2(Mxy,D2_plus-Myy);
end

