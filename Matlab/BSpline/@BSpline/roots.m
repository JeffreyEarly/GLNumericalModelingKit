function values = roots(spline)
%ROOTS Roots of BSpline for its domain

values = [];
scale = factorial((spline.K-1):-1:0);
C = spline.x_std*spline.C;
C(:,end) = C(:,end) + spline.x_mean;
t_pp = spline.t_pp;

for iBin=1:size(spline.C,1)
    r = roots(C(iBin,:)./scale);
    [~,I] = find( r >= 0 | r<= (t_pp(iBin+1)-t_pp(iBin)) );
    if ~isempty(I)
        values = cat(1,values,r(I)+t_pp(iBin));
    end
end

values = sort(values);