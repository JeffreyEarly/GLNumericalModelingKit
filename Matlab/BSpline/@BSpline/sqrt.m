function splinesqrt = sqrt(spline,constraints)
%SQRT Power of a BSpline
if exist('constraints','var')
    splinesqrt = power(spline,1/2,constraints);
else
    splinesqrt = power(spline,1/2);
end

