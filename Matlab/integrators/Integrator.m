classdef Integrator < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        stepSize
        currentTime
        totalIterations
        
        fFromTYVector
        currentY
        
        nReps
        nDims
        
        F
    end
    
    methods
        function self = Integrator( f, y0, dt )
            try
                f0 = feval(f,0,y0);
            catch theError
                msg = ['Unable to evaluate the ODEFUN at t0,y0. ',theError];
                error(msg);
            end
            
            [self.nReps, self.nDims] = size(y0);
            
            if ~isequal(size(y0),size(f0))
                error('Inconsistent sizes of Y0 and f(t0,y0).');
            end
            
            self.currentTime = 0;
            self.stepSize = dt;
            self.totalIterations = 0;
            self.fFromTYVector = f;
            self.currentY = y0;
            
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
                self.currentY = self.StepForward(self.currentY,self.currentTime,self.stepSize);
                self.currentTime = self.currentTime + self.stepSize;
                self.totalIterations = self.totalIterations + 1;
            end
            
            y = self.currentY;
        end
        
        function yo = StepForward(self,yi,t,dt)
            self.F(:,:,1) = feval(self.fFromTYVector,t,yi);
            self.F(:,:,2) = feval(self.fFromTYVector,t+0.5*dt,yi+0.5*dt*self.F(:,:,1));
            self.F(:,:,3) = feval(self.fFromTYVector,t+0.5*dt,yi+0.5*dt*self.F(:,:,2));
            self.F(:,:,4) = feval(self.fFromTYVector,t+dt,yi+dt*self.F(:,:,3));
            
            yo = yi + (dt/6)*(self.F(:,:,1) + 2*self.F(:,:,2) + 2*self.F(:,:,3) + self.F(:,:,4));
        end
    end
    
end

