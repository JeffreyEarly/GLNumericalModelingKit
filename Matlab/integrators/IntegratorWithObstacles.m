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
            
            if any( overlaps(obstacles)-eye(length(obstacles)))
                error('You have overlapping polygons! This is not allowed. To proceed use union to combine them.');
            end
                       
            self@Integrator( f, y0, dt);
            
            self.ymin = ymin;
            self.ymax = ymax;
            self.kappa = kappa;
            self.isPeriodic = isPeriodic;

            self.obstacles = obstacles;
            bufferRadius = 7*max(sqrt(2*kappa*dt));
            fprintf('based on your step size, we have to set the obstacle buffer radius to %f\n',bufferRadius);
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
    
            particlesToCheck = (1:size(yo,1))';
            nLoops = 0;
            while ~isempty(particlesToCheck)
                nLoops = nLoops + 1;

                % wrap the output for the particles we need to check.
                for i=1:self.nDims
                    yo(particlesToCheck,i) = yi(particlesToCheck,i) + dy(particlesToCheck,i);
                    if self.isPeriodic(i) == 1
                        yo(particlesToCheck,i) = mod(yo(particlesToCheck,i)-self.ymin(i),self.ymax(i)-self.ymin(i)) + self.ymin(i);
                    end
                end
                
                didReflectOnSomeObject = false(size(particlesToCheck));
                for iObstacle = 1:length(self.obstacles)
%                     if nLoops > 9
%                         fprintf('bail\n');
%                         figure, plot(self.obstacles(iObstacle))
%                         hold on, scatter(yi(particlesToCheck,1),yi(particlesToCheck,2),'filled')
%                         hold on, scatter(yo(particlesToCheck,1),yo(particlesToCheck,2),'filled')
%                     end
                    
                    [yi_new,dy_new,didReflect] = self.UpdateWithReflection(self.obstacles(iObstacle),yi(particlesToCheck,:),yo(particlesToCheck,:),dy(particlesToCheck,:));
                    
%                     if nLoops > 9
%                         hold on, scatter(yi_new(:,1)+dy_new(:,1),yi_new(:,2)+dy_new(:,2),'filled')
%                     end
                    
                    if any(didReflect)
                        reflectedIndices = particlesToCheck(didReflect);
                        yi(reflectedIndices,:) = yi_new(didReflect,:);
                        dy(reflectedIndices,:) = dy_new(didReflect,:);
                        didReflectOnSomeObject = didReflectOnSomeObject | didReflect;
                    end
                end
                particlesToCheck(~didReflectOnSomeObject) = [];
                
                if nLoops > 20
                    warning('Too many reflection loops. Something bad must of happened.')
                    particlesToCheck = [];
                end
            end
            
            yo = yi + dy;
        end
        

        function [yi,dy,didReflect] = UpdateWithReflection(self,obstacle,yi,yo_wrapped,dy)
            % yi and dy will contain updated values for indices where
            % didReflect == true.
            
            % defined using the increment, yi will be wrapped
            % the same number of times as yo.
            yi_wrapped = yo_wrapped - dy;
            [yo_fixed,p_intersect,didReflect] = IntegratorWithObstacles.Reflect(obstacle,yi_wrapped,yo_wrapped);
            if any(didReflect)
                dintersect = p_intersect(didReflect,:)- yi_wrapped(didReflect,:);
                dreflection = yo_fixed(didReflect,:) - p_intersect(didReflect,:);
                
                % The 'initial' point needs to becomes the
                % point of reflection
                yi(didReflect,:) = yi(didReflect,:) + dintersect;
                dy(didReflect,:) = dreflection;
            end
        end
        
    end
    
    methods (Static)
        function [p_f_new,p_intersection,didReflect] = Reflect(polygon,p_i,p_f)
            % Given line segment p_i-p_f, reflect it off the polygon if it
            % tries to cross inside. This only reflects *once*, so you need
            % to call this function until it stops reflecting.
            %
            % polygon is a polshape object
            % p_i is the initial position
            % p_f is the (proposed) final position
            %
            % What are the conditions for reflection?
            %
            % Define a "crossing" as p_i-p_f intersects a vertex. Note that
            % p_i or p_f being *on* a vertex counts as well.
            %
            % Easy cases first:
            % 1. p_i outside, p_f inside, 1..3..5.. all crossings have dist > 0 (reflect at first crossing).
            % 2. p_i outside, p_f outside, 2..4..6.. all crossings have dist > 0 (reflect at first crossing)
            % Tricky cases now:
            % 3. p_i on border, p_f inside, ODD crossings (reflect at first crossing, particle went the wrong way at the start)
            % 4. p_i on border, p_f inside, EVEN crossings (reflect at *second* crossing, particle was out, but came back in). *** if p_f on border, it's a NO-OP
            % 5. p_i on border, p_f outside, 1 cross --- NO-OP, particle  stayed ouside
            % 6. p_i on border, p_f outside, ODD (>1) crossings (reflect at *second* crossing, particle was out, but came back in)
            % 7. p_i on border, p_f outside, EVEN crossings (reflect at first crossing, particle went the wrong way at the start)
            % 
            % Possible outcomes:
            % 1. NO-OP
            % 2. Reflect at first crossing
            % 3. Reflect at second crossing. This case only occurs if p_i
            % is on a vertex. So we only ever need to store the first
            % non-zero crossing.
            
            % output
            didReflect = false(size(p_i,1),1);
            p_f_new = p_f;
            p_intersection = zeros(size(p_f));
            
            % Information needed to determine whether or not we reflect
            isPfInterior = isinterior(polygon,p_f(:,1),p_f(:,2));
            isPiOnBorder = false(size(p_i,1),1);
            numCrossings = zeros(size(p_i,1),1);
            
            % Information needed to do the actual reflection
            first_non_zero_crossing = zeros(size(p_f));
            dist_to_first_nonzero_crossing = nan(size(p_i,1),1);
            
            zeroThreshold = 0; % this should be a relative value.
            for i=1:length(polygon.Vertices)
                if i == length(polygon.Vertices)
                    polyedge = [polygon.Vertices(end,1:2); polygon.Vertices(1,1:2)];
                else
                    polyedge = polygon.Vertices(i:(i+1),1:2);
                end
                [doesIntersect,p_intersect,r2] = IntegratorWithObstacles.DoesIntersect( p_i, p_f, polyedge(1,:), polyedge(2,:) );
                
                numCrossings = numCrossings + doesIntersect;
                isPiOnBorder = isPiOnBorder | r2 <= zeroThreshold;
                
                % *** New intersection to record if ***
                % it intersects, isn't on the border, and we don't previously have a point recorded, OR
                % it isn't on the border and is closer
                newIntersection = (doesIntersect & r2 > zeroThreshold & isnan(dist_to_first_nonzero_crossing)) | (r2 > zeroThreshold & r2 < dist_to_first_nonzero_crossing);
                if any(newIntersection)
                   dist_to_first_nonzero_crossing(newIntersection) = r2(newIntersection);
                   first_non_zero_crossing(newIntersection,:) = p_intersect(newIntersection,:);
                end
                
                % Update new possible reflection.
                % Reflect if we crossed a boundary, except case 5 above
                shouldReflect = doesIntersect & numCrossings > 0 & ~(isPiOnBorder & ~isPfInterior & numCrossings == 1);
                didReflect = didReflect | shouldReflect;
                if any(shouldReflect)
                    p_intersection(doesIntersect,:) = first_non_zero_crossing(doesIntersect,:);
                    
                    % Reflection point is p_i in cases 3 and 7
                    piIsReflectionPoint = shouldReflect & ( (isPiOnBorder & isPfInterior & mod(numCrossings,2) == 1) | (isPiOnBorder & ~isPfInterior & mod(numCrossings,2) == 0) );
                    p_intersection(piIsReflectionPoint,:) = p_i(piIsReflectionPoint,:);
                    
                    % treat the point of intersection as the origin. This means the
                    % interior point is
                    px = p_f(shouldReflect,1) - p_intersect(shouldReflect,1);
                    py = p_f(shouldReflect,2) - p_intersect(shouldReflect,2);
                    
                    dx = polyedge(2,1)-polyedge(1,1);
                    dy = polyedge(2,2)-polyedge(1,2);
                    
                    a = dx*dx-dy*dy;
                    b = 2*dx*dy;
                    d = dx*dx+dy*dy;
                    
                    p_f_new(shouldReflect,1) = (a*px + b*py)/d + p_intersect(shouldReflect,1);
                    p_f_new(shouldReflect,2) = (b*px - a*py)/d + p_intersect(shouldReflect,2);
                end
            end
%             if any(isPfInterior & ~didReflect)
%                fprintf('What happened here?'); 
%             end
        end
        
        function [doesIntersect,p_intersect,r2] = DoesIntersect(p_i,p_f,vi,vf)
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
            r2 = nan(size(p_i,1),1);
            doesIntersect = false(size(p_i,1),1);
            t = zeros(size(p_i,1),1);
            u = zeros(size(p_i,1),1);
            
            valid = denom ~= 0;
            t(valid) = a(valid)./denom(valid);
            u(valid) = -b(valid)./denom(valid);
            
            % if t==0, then the initial point is *on* a vertex. This should
            % not be used for a reflection
            doesIntersect(valid) = logical(t(valid) >= 0 & t(valid) <= 1 & u(valid) >= 0 & u(valid) <= 1);
            if any(doesIntersect)
                dx = - 0.999*t(doesIntersect).*x1x2(doesIntersect);
                dy = - 0.999*t(doesIntersect).*y1y2(doesIntersect);
                r2(doesIntersect) = dx.^2+dy.^2;
                p_intersect(doesIntersect,1) = p_i(doesIntersect,1) + dx;
                p_intersect(doesIntersect,2) = p_i(doesIntersect,2) + dy;
            end
        end
    end
end

