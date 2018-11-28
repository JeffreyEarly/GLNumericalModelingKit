function flag = IsEvenlySampled(t)
% Checks the sampling rate of t.
% Returns 1 if the data is evenly sampled (a single unique dt)
% Returns 2 if the data is sampled with multiples of a unique dt
% Return 0 otherwise
unique_dt = unique(diff(t));
if length(unique_dt) == 1
    flag = 1;
else
    dt_multiples = unique_dt/min(unique_dt);
    if all( dt_multiples-1.0 < 1e-7 )
        flag = 1;
    elseif all(mod(dt_multiples,1.0) < 0.01)
        flag = 2;
    else
        flag = 0;
    end
end
end