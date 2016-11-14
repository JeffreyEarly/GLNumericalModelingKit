function [smoothVec] = RunningAverage( vec, smoothness )
%	Returns a running average of vec. Tries to handle end points well.

smoothVec = RunningFilter(vec,smoothness,'mean');



