classdef IntegratorEulerMaruyama < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        stepSize
        currentTime
        totalIterations
        
        fFromTYVector
        gFromTYVector
        currentY
        lastY
        lastTime;
        
        nReps
        nDims
        
        F
    end
    
    methods
        function self = IntegratorEulerMaruyama( f, g, y0, dt )
            try
                f0 = feval(f,0,y0);
            catch theError
                msg = ['Unable to evaluate f(t0,y0). ',theError];
                error(msg);
            end
            
            try
                g0 = feval(g,0,y0);
            catch theError
                msg = ['Unable to evaluate g(t0,y0). ',theError];
                error(msg);
            end
            
            [self.nReps, self.nDims] = size(y0);
            
            if ~isequal(size(y0),size(f0),size(g0))
                error('Inconsistent sizes of y0 and f(t0,y0) and g(t0,y0).');
            end
            
            self.currentTime = 0;
            self.stepSize = dt;
            self.totalIterations = 0;
            self.fFromTYVector = f;
            self.gFromTYVector = g;
            self.currentY = y0;
            self.lastY = [];
            self.lastTime = [];
            
            self.F = zeros(self.nReps,self.nDims,4);
        end
        
        function [y, t] = IntegrateAlongDimension(self,time)
            y = zeros(self.nReps,self.nDims,length(time));
            t = zeros(size(time));
            for i=1:length(time)
                p = self.StepForwardToTime(time(i));
                y(:,:,i) = p;
                t(i) = self.currentTime;
            end
        end
        
        function y = StepForwardToTime(self, time )
            while self.currentTime < time
                self.lastY = self.currentY;
                self.lastTime = self.currentTime;
                self.currentY = self.StepForward(self.currentY,self.currentTime,self.stepSize);
                self.currentTime = self.currentTime + self.stepSize;
                self.totalIterations = self.totalIterations + 1;
            end
            
            if isempty(self.lastY)
                y = self.currentY;
            else
                % Euler-Maruyama is 1st order, which means we can do linear
                % interpolation without formally losing accuracy.
                alpha = (self.currentTime-time)/(self.currentTime-self.lastTime);
                y = alpha*self.lastY + (1-alpha)*self.currentY;
            end
        end
        
        function yo = StepForward(self,yi,t,dt)
            yo = yi + dt*feval(self.fFromTYVector,t,yi) + sqrt(dt)*randn(size(yi)).*feval(self.gFromTYVector,t,yi);
        end
    end
    
end

