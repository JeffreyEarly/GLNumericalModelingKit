function [p_f_new,p_intersection,didReflect] = Reflect(polygon,p_i,p_f)
            % Ugly, optimized version
            %
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
            isPfInterior = IntegratorWithObstacles.isInterior( polygon,p_f(:,1),p_f(:,2));
            isPiOnBorder = false(size(p_i,1),1);
            numCrossings = zeros(size(p_i,1),1);
                        
            % Information needed to do the actual reflection
            first_non_zero_crossing = zeros(size(p_f));
            dist_to_first_nonzero_crossing = nan(size(p_i,1),1);
            
            zeroThreshold = 0; % this should be a relative value.
            doesIntersect = false(size(p_i,1),1);
            p_intersect = nan(size(p_i));
            r2 = nan(size(p_i,1),1);
            for i=1:(length(polygon.Vertices)-1)
%                 [doesIntersect,p_intersect,r2] = IntegratorWithObstacles.DoesIntersect( p_i, p_f, polygon.Vertices(i,1:2), polygon.Vertices(i+1,1:2), doesIntersect,p_intersect,r2 );
                                
                x1x2 = p_i(:,1)-p_f(:,1);   % x1 - x2
                y3y4 = -polygon.dy(i);      % y3 - y4
                y1y2 = p_i(:,2)-p_f(:,2);   % y1 - y2
                x3x4 = -polygon.dx(i);      % x3 - x4
                x1x3 = p_i(:,1)-polygon.Vertices(i,1);    % x1 - x3
                y1y3 = p_i(:,2)-polygon.Vertices(i,2);    % y1 - y3
                
                denom = x1x2 * y3y4 - y1y2 * x3x4;
                a = x1x3 .* y3y4 - y1y3 .* x3x4;
                b = x1x2 .* y1y3 - y1y2 .* x1x3;
                
                t = a./denom;
                u = -b./denom;
                
                % if t==0, then the initial point is *on* a vertex. This should
                % not be used for a reflection
                doesIntersect = t >= 0 & t <= 1 & u >= 0 & u <= 1 & denom ~= 0;
                if any(doesIntersect)
                    dx = - 0.999*t(doesIntersect).*x1x2(doesIntersect);
                    dy = - 0.999*t(doesIntersect).*y1y2(doesIntersect);
                    r2(doesIntersect) = dx.^2+dy.^2;
                    p_intersect(doesIntersect,1) = p_i(doesIntersect,1) + dx;
                    p_intersect(doesIntersect,2) = p_i(doesIntersect,2) + dy;
                end
                
                
                
                
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
                    
                    dx = polygon.Vertices(i+1,1)-polygon.Vertices(i,1);
                    dy = polygon.Vertices(i+1,2)-polygon.Vertices(i,2);
                    
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