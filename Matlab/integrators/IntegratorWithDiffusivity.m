classdef IntegratorWithDiffusivity < Integrator
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        ymin
        ymax
        kappa
        
        diffusivityFlux
    end
    
    methods
        function self = IntegratorWithDiffusivity( f, y0, dt, kappa, ymin, ymax )
            [~, nDims] = size(y0);
            
            if length(kappa) == 1 && nDims > 1 % if it's a scalar, make it a vector
                kappa = kappa*ones(1,nDims);
            end
            if ~isequal(size(kappa,2),nDims)
                error('Inconsistent sizes of Y0 and kappa.');
            end
            
            if ~isequal([1 nDims],size(ymin)) || ~isequal([1 nDims],size(ymax))
                error('Inconsistent sizes of Y0 and ymin/ymax.');
            end
                       
            self@Integrator( f, y0, dt);
            
            self.ymin = ymin;
            self.ymax = ymax;
            self.kappa = kappa;
            
            self.diffusivityFlux = cell(length(ymin),1);
            for i = 1:length(ymin)
                if kappa(i) == 0
                    % no-op
                    self.diffusivityFlux{i} = @(x,dt) zeros(size(x));
                elseif ymin(i) == -Inf && isfinite(ymax(i))
                    % bounded above
                    b = ymax(i);
                    s = sqrt(2*kappa(i));
                    self.diffusivityFlux{i} = @(x,dt) -abs((s/sqrt(dt))*randn(size(x)) - (b-x)/dt) + (b-x)/dt;
                elseif isfinite(ymin(i)) && ymax(i) == Inf
                    % bounded below
                    a = ymin(i);
                    s = sqrt(2*kappa(i));
                    self.diffusivityFlux{i} = @(x,dt) abs((s/sqrt(dt))*randn(size(x)) - (a-x)/dt) + (a-x)/dt;
                elseif isfinite(ymin(i)) && isfinite(ymax(i))
                    % bounded above and below
                    a = ymin(i);
                    b = ymax(i);
                    s = sqrt(2*kappa(i));
                    self.diffusivityFlux{i} = @(x,dt) abs( mod( (s/sqrt(dt))*randn(size(x))-(b-x)/dt, 2*(b-a)/dt) - (b-a)/dt ) + (a-x)/dt;
                elseif ymin(i) == -Inf && ymax(i) == Inf
                    % unbounded
                    self.diffusivityFlux{i} = @(x,dt) sqrt(2*kappa(i)/dt)*randn(size(x));
                else
                    error('Invalid diffusion boundary')
                end
            end
        end
        
        
        function yo = StepForward(self,yi,t,dt)            
            yo = StepForward@Integrator(self,yi,t,dt);
            for i=1:self.nDims
               yo(:,i) = yo(:,i) + dt*self.diffusivityFlux{i}(yo(:,i),dt); 
            end
        end
    end
    
end

