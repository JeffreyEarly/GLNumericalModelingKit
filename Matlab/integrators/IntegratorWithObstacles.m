classdef IntegratorWithObstacles < Integrator
    %UNTITLED 2D only
    
    properties
        ymin
        ymax
        kappa
        isPeriodic = [0 0];
        
        obstacles
        obstacleBuffers
        obstacleBoundingBoxes
        
        diffusivityFlux
    end
    
    methods
        function self = IntegratorWithObstacles( f, y0, dt, kappa, ymin, ymax, obstacles, isPeriodic )
            [~, nDims] = size(y0);
            
            if length(kappa) == 1 && nDims > 1 % if it's a scalar, make it a vector
                kappa = kappa*ones(1,nDims);
            end
            if ~isequal(size(kappa,2),nDims)
                error('Inconsistent sizes of Y0 and kappa.');
            end
            if ~isequal(length(isPeriodic),nDims)
                error('Inconsistent sizes of isPeriodic and Y0.');
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
            self.isPeriodic = isPeriodic;

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
                elseif ymin(i) == -Inf && isfinite(ymax(i)) && isPeriodic(i) == 0
                    % bounded above
                    b = ymax(i);
                    s = sqrt(2*kappa(i));
                    self.diffusivityFlux{i} = @(x,dt) -abs((s/sqrt(dt))*randn(size(x)) - (b-x)/dt) + (b-x)/dt;
                elseif isfinite(ymin(i)) && ymax(i) == Inf && isPeriodic(i) == 0
                    % bounded below
                    a = ymin(i);
                    s = sqrt(2*kappa(i));
                    self.diffusivityFlux{i} = @(x,dt) abs((s/sqrt(dt))*randn(size(x)) - (a-x)/dt) + (a-x)/dt;
                elseif isfinite(ymin(i)) && isfinite(ymax(i)) && isPeriodic(i) == 0
                    % bounded above and below
                    a = ymin(i);
                    b = ymax(i);
                    s = sqrt(2*kappa(i));
                    self.diffusivityFlux{i} = @(x,dt) abs( mod( (s/sqrt(dt))*randn(size(x))-(b-x)/dt, 2*(b-a)/dt) - (b-a)/dt ) + (a-x)/dt;
                else
                    % unbounded
                    self.diffusivityFlux{i} = @(x,dt) sqrt(2*kappa(i)/dt)*randn(size(x));
%                 else
%                     error('Invalid diffusion boundary')
                end
            end
        end
        
        function yo = StepForward(self,yi,t,dt)            
            yo = StepForward@Integrator(self,yi,t,dt);
            
            % Obstacles get tricky when dealing with periodic boundaries.
            % Conceptually, the way we've decided to tackle this is to
            % check if the wrapped final position is in an obstacle,
            % reflect, and then compute the increment from the unwrapped
            % value.
            dy = yo-yi; % increment from RK4
            for i=1:self.nDims % plus the increment from diffusivity
                dy(:,i) = dy(:,i) + dt*self.diffusivityFlux{i}(yo(:,i),dt);
            end
            
            
            
            % Let's check all particles to see if the end up in an
            % obstacle.
            particlesToCheck = 1:size(yo,1);
            while ~isempty(particlesToCheck)
                % wrap the output for the particles we need to check.
                for i=1:self.nDims
                    yo(particlesToCheck,i) = yi(particlesToCheck,i) + dy(particlesToCheck,i);
                    if self.isPeriodic(i) == 1
                        yo(particlesToCheck,i) = mod(yo(particlesToCheck,i)-self.ymin(i),self.ymax(i)-self.ymin(i)) + self.ymin(i);
                    end
                end
                for iObstacle = 1:length(self.obstacles)
                    isInside = isinterior(self.obstacles(iObstacle),yo(particlesToCheck,1),yo(particlesToCheck,2));
                    particlesToCheck(~isInside) = [];
                    if ~isempty(particlesToCheck)
                        % defined using the increment, yi will be wrapped
                        % the same number of times as yo.
                        yi_wrapped = yo(particlesToCheck,:) - dy(particlesToCheck,:);
                        [yo_fixed,p_intersect,didReflect] = IntegratorWithObstacles.Reflect(self.obstacles(iObstacle),yi_wrapped,yo(particlesToCheck,:));
                        if any(didReflect)
                            dintersect = p_intersect(didReflect,:)- yi_wrapped; 
                            dreflection = yo_fixed(didReflect,:) - p_intersect(didReflect,:);
                            
                            % The 'initial' point needs to becomes the
                            % point of reflection
                            yi(particlesToCheck(didReflect),:) = yi(particlesToCheck(didReflect),:) + dintersect;
                            dy(particlesToCheck(didReflect),:) = dreflection;
                        end
                    end
                end
            end
            
            yo = yi + dy;
        end
        

        
        
    end
    
    methods (Static)
        function [p_f,p_intersection,didReflect] = Reflect(polygon,p_i,p_f)
            % polgon is a polshape object
            % p_i is the initial position
            % p_f is the (proposed) final position
            
            didReflect = false(size(p_i,1),1);
            p_intersection = zeros(size(p_f));
            for i=1:length(polygon.Vertices)
                if i == length(polygon.Vertices)
                    polyedge = [polygon.Vertices(end,1:2); polygon.Vertices(1,1:2)];
                else
                    polyedge = polygon.Vertices(i:(i+1),1:2);
                end
                [doesIntersect,p_intersect] = IntegratorWithObstacles.DoesIntersect( p_i, p_f, polyedge(1,:), polyedge(2,:) );
                
                didReflect = didReflect | doesIntersect;
                if any(doesIntersect)
                    % We should do a distance check! Feflect on the closest
                    % point only.
                    p_intersection(doesIntersect,:) = p_intersect(doesIntersect,:);
                    % treat the point of intersection as the origin. This means the
                    % interior point is
                    px = p_f(doesIntersect,1) - p_intersect(doesIntersect,1);
                    py = p_f(doesIntersect,2) - p_intersect(doesIntersect,2);
                    
                    dx = polyedge(2,1)-polyedge(1,1);
                    dy = polyedge(2,2)-polyedge(1,2);
                    
                    a = dx*dx-dy*dy;
                    b = 2*dx*dy;
                    d = dx*dx+dy*dy;
                    
                    p_f(doesIntersect,1) = (a*px + b*py)/d + p_intersect(doesIntersect,1);
                    p_f(doesIntersect,2) = (b*px - a*py)/d + p_intersect(doesIntersect,2);
                end
            end
        end
        
        function [doesIntersect,p_intersect] = DoesIntersect(p_i,p_f,vi,vf)
            % somewhat vectorized intersection algorithm
            % pi and pf are the initial and final positions, given as [x y]
            % with as many rows as you want.
            % vi and vf are the initial and final vertices of the polygon,
            % in the same format, but only one row.
            % https://en.wikipedia.org/wiki/Line–line_intersection
            
            % flag will indicate whether the line intersects the vertex.
            
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
            doesIntersect = false(size(p_i,1),1);
            t = zeros(size(p_i,1),1);
            u = zeros(size(p_i,1),1);
            
            valid = denom ~= 0;
            t(valid) = a(valid)./denom(valid);
            u(valid) = -b(valid)./denom(valid);
            
            doesIntersect(valid) = logical(t(valid) >= 0 & t(valid) <= 1 & u(valid) >= 0 & u(valid) <= 1);
            if any(doesIntersect)
                p_intersect(doesIntersect,1) = p_i(doesIntersect,1) - t(doesIntersect).*x1x2(doesIntersect);
                p_intersect(doesIntersect,2) = p_i(doesIntersect,2) - t(doesIntersect).*y1y2(doesIntersect);
            end
        end
    end
end

