classdef IntegratorWithObstacles < Integrator
    %UNTITLED 2D only
    
    properties
        ymin
        ymax
        kappa
        
        obstacles
        obstacleBuffers
        obstacleBoundingBoxes
        
        diffusivityFlux
    end
    
    methods
        function self = IntegratorWithObstacles( f, y0, dt, kappa, ymin, ymax, obstacles )
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
            
            if nDims ~= 2
                error('This integrator is only valid for 2D.')
            end
                       
            self@Integrator( f, y0, dt);
            
            self.ymin = ymin;
            self.ymax = ymax;
            self.kappa = kappa;

            self.obstacles = obstacles;
            bufferRadius = 5*max(sqrt(2*kappa*dt));
            if ~isempty(self.obstacles)
                self.obstacleBuffers = polybuffer(self.obstacles,bufferRadius);
            end
            
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

            for iObstacle = 1:length(self.obstacles)
                isInside = isinterior(self.obstacleBuffers(iObstacle),yo(:,1),yo(:,2));
                if any(isInside)
                    [yo(isInside,:),numReflections] = IntegratorWithObstacles.Reflect(self.obstacles(iObstacle),yi(isInside,:),yo(isInside,:));
                    fprintf('Found %d particles inside our buffer; reflected %d of them.\n',sum(isInside),numReflections);
                end
            end
        end
        

        
        
    end
    
    methods (Static)
        function [p_f,numReflections] = Reflect(polgon,p_i,p_f)
            % polgon is a polshape object
            % p_i is the initial position
            % p_f is the (proposed) final position
            
            numReflections = 0;
            for i=1:length(polgon.Vertices)
                if i == length(polgon.Vertices)
                    polyedge = [polgon.Vertices(end,1:2); polgon.Vertices(1,1:2)];
                else
                    polyedge = polgon.Vertices(i:(i+1),1:2);
                end
                [flag,p_intersect] = IntegratorWithObstacles.DoesIntersect( p_i, p_f, polyedge(1,:), polyedge(2,:) );
                
                numReflections = numReflections + sum(flag);
                if any(flag)
                    % treat the point of intersection as the origin. This means the
                    % interior point is
                    px = p_f(flag,1) - p_intersect(flag,1);
                    py = p_f(flag,2) - p_intersect(flag,2);
                    
                    dx = polyedge(2,1)-polyedge(1,1);
                    dy = polyedge(2,2)-polyedge(1,2);
                    
                    a = dx*dx-dy*dy;
                    b = 2*dx*dy;
                    d = dx*dx+dy*dy;
                    
                    p_f(flag,1) = (a*px + b*py)/d + p_intersect(flag,1);
                    p_f(flag,2) = (b*px - a*py)/d + p_intersect(flag,2);
                end
            end
        end
        
        function [flag,p_intersect] = DoesIntersect(p_i,p_f,vi,vf)
            % somewhat vectorized intersection algorithm
            % pi and pf are the initial and final positions, given as [x y]
            % with as many rows as you want.
            % vi and vf are the initial and final vertices of the polygon,
            % in the same format, but only one row.
            % https://en.wikipedia.org/wiki/Line–line_intersection
            
            x1x2 = p_i(:,1)-p_f(:,1);   % x1 - x2
            y3y4 = vi(2)-vf(2);         % y3 - y4
            y1y2 = p_i(:,2)-p_f(:,2);   % y1 - y2
            x3x4 = vi(1)-vf(1);         % x3 - x4
            x1x3 = p_i(:,1)-vi(:,1);    % x1 - x3
            y1y3 = p_i(:,2)-vi(:,2);    % y1 - y3
            
            denom = x1x2 * y3y4 - y1y2 * x3x4;
            a = x1x3 .* y3y4 - y1y3 .* x3x4;
            b = x1x2 .* y1y3 - y1y2 .* x1x3;
            
            p_intersect = nan(size(p_i));
            flag = false(size(p_i,1),1);
            t = zeros(size(p_i,1),1);
            u = zeros(size(p_i,1),1);
            
            valid = denom ~= 0;
            t(valid) = a(valid)./denom(valid);
            u(valid) = -b(valid)./denom(valid);
            
            flag(valid) = logical(t >= 0 & t <= 1 & u >= 0 & u <= 1);
            if any(flag)
                p_intersect(flag,1) = p_i(flag,1) - t(flag).*x1x2(flag);
                p_intersect(flag,2) = p_i(flag,2) - t(flag).*y1y2(flag);
            end
        end
    end
end

