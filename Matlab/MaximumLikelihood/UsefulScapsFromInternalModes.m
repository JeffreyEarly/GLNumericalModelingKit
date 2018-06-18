       function rho_zCheb = FitGriddedDensityDataToMonotonicFunctionWithSplines(z,rho)
            % useful notes about condition number in least squares
            % http://www.math.uchicago.edu/~may/REU2012/REUPapers/Lee.pdf
            n = length(z);
            Lz = max(z)-min(z);
            T = InternalModesSpectral.ChebyshevPolynomialsOnGrid(z,n);
            maxPolys = InternalModes.NumberOfWellConditionedModes(T);
            rho_zCheb = T(:,1:maxPolys)\rho;
            
            % dense grid, upon which to evaluate derivatives
            n_dense = 10*maxPolys;
            
            rho_z_zCheb = InternalModesSpectral.DifferentiateChebyshevVector(rho_zCheb);
            rho_z_zCheb = cat(1,rho_z_zCheb,zeros(n_dense-length(rho_z_zCheb),1));
            rho_z = InternalModesSpectral.ifct(rho_z_zCheb);
            
            if any(rho_z < 0)
                fprintf('Fit gridded data to %d modes. Bouyancy frequency goes negative; attempting to constrain.\n',maxPolys)
                
                meanRho = mean(rho);
                rho = rho-meanRho;
                
                sigma_rho = sqrt(mean((rho-mean(rho)).^2));
                sigma_drhodz = sqrt(mean((diff(rho)./diff(z)).^2));
                lambda = n*sigma_rho*sigma_rho/(1e-8*Lz*sigma_drhodz*sigma_drhodz);
                
                % linear function, describing the likelihood of a derivative occurring.
                min_rho_z = -1e-4*sigma_drhodz;
                max_rho_z = 1e4*max(diff(rho)./diff(z));
                min_lambda = lambda;
                max_lambda = 1e4*min_lambda;
                slope = (max_lambda-min_lambda)/(max_rho_z-min_rho_z);
                intercept = min_lambda - slope*min_rho_z;
                
                w = @(rho_z) (slope*rho_z+intercept).*(rho_z > min_rho_z);
                
                z_dense = linspace(min(z),max(z),n_dense)';
                dz = z_dense(2)-z_dense(1);
                
                iPoly = 3;
                K = 5;
                
                z_lin = linspace(min(z),max(z),length(z))';
                t_knot = BSpline.KnotPointsForSplines(z,K,iPoly);
                T = BSpline.Spline( z, t_knot, K );
                B = BSpline.Spline( z_dense, t_knot, K, 1 );
                Tz = squeeze(B(:,:,2));
                
                f = @(m) ((rho-T*m).')*(rho-T*m) + ((Tz*m-min_rho_z).')*(diag(w(Tz*m).*dz))*(Tz*m-min_rho_z);
                
                % First we use the Newton iteration across some low-order
                % fits to find a good starting point;
                
                iPoly = 3;
                f = @(m) ((rho-T*m).')*(rho-T*m) + ((Tz*m-min_rho_z).')*(diag(w(Tz*m).*dz))*(Tz*m-min_rho_z);
                rho_zCheb = InternalModesSpectral.FindMinimum(T\rho,T,Tz,w,dz,rho,min_rho_z);
                
                rho_zChebMin = rho_zCheb;
                fmin0 = f(rho_zChebMin);
                k = iPoly;
                AICc_min = 2*k + n*log(fmin0) + (2*k*k+2*k)/(n-k-1);
                rho_min = BSpline(z,rho,K,t_knot,rho_zCheb);
                for iPoly = 4:10
                    t_knot = BSpline.KnotPointsForSplines(z,K,iPoly);
                    T = BSpline.Spline( z, t_knot, K );
                    B = BSpline.Spline( z_dense, t_knot, K, 1 );
                    Tz = squeeze(B(:,:,2));
                    f = @(m) ((rho-T*m).')*(rho-T*m) + ((Tz*m-min_rho_z).')*(diag(w(Tz*m).*dz))*(Tz*m-min_rho_z);
                    rho_zCheb = InternalModesSpectral.FindMinimum(T\rho,T,Tz,w,dz,rho,min_rho_z);
                    
                    [rho_zCheb, fnew] = fminsearch(f,rho_zCheb,optimset('TolX',1e-4));
                    %                     fnew = f(rho_zCheb);
                    
                    k = size(T,2);
                    AIC = 2*k + n*log(fnew);
                    AICc = AIC + (2*k*k+2*k)/(n-k-1);
                    fprintf('iPoly: %d, fmin: %f, AIC: %f, AICc: %f\n', k, fnew, AIC, AICc);
                    
                    if AICc < AICc_min
                        fmin0 = fnew;
                        AICc_min = AICc;
                        rho_zChebMin = rho_zCheb;
                        rho_min = BSpline(z,rho,K,t_knot,rho_zCheb);
                    end
                end
                fprintf('Initial iteration found %d Polys: fmin = %f \n', length(rho_zChebMin), fmin0);
                    
                f = @(x) rho_min(x) + meanRho;
                
                [self.zLobatto, rho_zCheb] = InternalModesSpectral.ProjectOntoChebyshevPolynomialsWithTolerance([min(z) max(z)], f, 1e-16);
                
            else
                fprintf('Fit gridded data to %d modes. Bouyancy frequency stays positive.\n',maxPolys)
                return
            end
            
            %             self.zLobatto = (Lz/2)*( cos(((0:n-1)')*pi/(n-1)) + 1) + min(z);
            %             self.rho_zCheb = reshape(T\rho,[],1);
            %             self.rho_zCheb = cat(1,self.rho_zCheb,zeros(n-length(self.rho_zCheb),1));
            %             self.rho_zLobatto = InternalModesSpectral.ifct(self.rho_zCheb);
        end
        
        function rho_zCheb = FitGriddedDensityDataToMonotonicFunction(z,rho)
            % useful notes about condition number in least squares
            % http://www.math.uchicago.edu/~may/REU2012/REUPapers/Lee.pdf
            n = length(z);
            Lz = max(z)-min(z);
            T = InternalModesSpectral.ChebyshevPolynomialsOnGrid(z,n);
            maxPolys = InternalModes.NumberOfWellConditionedModes(T);
            rho_zCheb = T(:,1:maxPolys)\rho;
            
            % dense grid, upon which to evaluate derivatives
            n_dense = 10*maxPolys;
            
            rho_z_zCheb = InternalModesSpectral.DifferentiateChebyshevVector(rho_zCheb);
            rho_z_zCheb = cat(1,rho_z_zCheb,zeros(n_dense-length(rho_z_zCheb),1));
            rho_z = InternalModesSpectral.ifct(rho_z_zCheb);
            
            if any(rho_z < 0) 
                fprintf('Fit gridded data to %d modes. Bouyancy frequency goes negative; attempting to constrain.\n',maxPolys)
                
                sigma_rho = sqrt(mean((rho-mean(rho)).^2));
                sigma_drhodz = sqrt(mean((diff(rho)./diff(z)).^2));
                lambda = n*sigma_rho*sigma_rho/(1e-8*Lz*sigma_drhodz*sigma_drhodz);
                
                % linear function, describing the likelihood of a derivative occurring.
                min_rho_z = -1e-4*sigma_drhodz;
                max_rho_z = 1e4*max(diff(rho)./diff(z));
                min_lambda = lambda;
                max_lambda = 1e4*min_lambda;
                slope = (max_lambda-min_lambda)/(max_rho_z-min_rho_z);
                intercept = min_lambda - slope*min_rho_z;
                
                w = @(rho_z) (slope*rho_z+intercept).*(rho_z > min_rho_z);
                
                
                z_dense = ((max(z)-min(z))/2)*( cos(((0:n_dense-1)')*pi/(n_dense-1)) + 1) + min(z)';
                [~,Tz] = InternalModesSpectral.ChebyshevPolynomialsOnGrid(z_dense,n);
                Tz = Tz(:,1:maxPolys);
                
                % very rough integration weights
                dz = diff(z_dense);
                dz = abs(cat(1,dz(1)/2,(dz(1:end-1)+dz(2:end))/2, dz(1)/2));
                
                f = @(m) ((rho-T*m).')*(rho-T*m) + ((Tz*m-min_rho_z).')*(diag(w(Tz*m).*dz))*(Tz*m-min_rho_z);
                
                % First we use the Newton iteration across some low-order
                % fits to find a good starting point;
                
                iPoly = 3;
                f = @(m) ((rho-T(:,1:iPoly)*m).')*(rho-T(:,1:iPoly)*m) + ((Tz(:,1:iPoly)*m-min_rho_z).')*(diag(w(Tz(:,1:iPoly)*m).*dz))*(Tz(:,1:iPoly)*m-min_rho_z);
                rho_zCheb = InternalModesSpectral.FindMinimum(T(:,1:iPoly)\rho,T(:,1:iPoly),Tz(:,1:iPoly),w,dz,rho,min_rho_z);
                
                rho_zChebMin = rho_zCheb;
                fmin0 = f(rho_zChebMin);
                k = iPoly;
                AICc_min = 2*k + n*log(fmin0) + (2*k*k+2*k)/(n-k-1);
                for iPoly = 4:10
                    f = @(m) ((rho-T(:,1:iPoly)*m).')*(rho-T(:,1:iPoly)*m) + ((Tz(:,1:iPoly)*m-min_rho_z).')*(diag(w(Tz(:,1:iPoly)*m).*dz))*(Tz(:,1:iPoly)*m-min_rho_z);
                    rho_zCheb = InternalModesSpectral.FindMinimum(T(:,1:iPoly)\rho,T(:,1:iPoly),Tz(:,1:iPoly),w,dz,rho,min_rho_z);
                    [rho_zCheb, fnew] = fminsearch(f,rho_zCheb,optimset('TolX',1e-4));
%                     fnew = f(rho_zCheb);

                    k = iPoly;
                    AIC = 2*k + n*log(fnew);
                    AICc = AIC + (2*k*k+2*k)/(n-k-1);
                    fprintf('iPoly: %d, fmin: %f, AIC: %f, AICc: %f\n', k, fnew, AIC, AICc);
                    
                    if AICc < AICc_min
                        fmin0 = fnew;
                        AICc_min = AICc;
                        rho_zChebMin = rho_zCheb;
                    end
                end
                fprintf('Initial iteration found %d Polys: fmin = %f \n', length(rho_zChebMin), fmin0);
                
                initialPolys = length(rho_zChebMin);
                dPoly = 4;
                rho_zCheb = rho_zChebMin;
%                 for iPoly = initialPolys:dPoly:(maxPolys/2)
%                     rho_zCheb = cat(1,rho_zCheb,zeros(iPoly-length(rho_zCheb),1));
%                     f = @(m) ((rho-T(:,1:iPoly)*m).')*(rho-T(:,1:iPoly)*m) + ((Tz(:,1:iPoly)*m-min_rho_z).')*(diag(w(Tz(:,1:iPoly)*m).*dz))*(Tz(:,1:iPoly)*m-min_rho_z);
%                     [rho_zCheb, fmin_final] = fminsearch(f,rho_zCheb,optimset('TolX',1e-4));
%                     fprintf('%d Polys: fmin = %f\n', iPoly, fmin_final);
%                 end
                
            else
                fprintf('Fit gridded data to %d modes. Bouyancy frequency stays positive.\n',maxPolys)
                return
            end
            
%             self.zLobatto = (Lz/2)*( cos(((0:n-1)')*pi/(n-1)) + 1) + min(z);
%             self.rho_zCheb = reshape(T\rho,[],1);
%             self.rho_zCheb = cat(1,self.rho_zCheb,zeros(n-length(self.rho_zCheb),1));
%             self.rho_zLobatto = InternalModesSpectral.ifct(self.rho_zCheb);
        end
        
        function m = FindMinimum(m0, T, Tz, w, dz, rho, min_rho_z)
            m = m0;
            TT = (T.')*T;
            
            max_dm = 1;
            iterations = 0;
            while max_dm > 1e-4 && iterations < 1000   
                rho_z = Tz*m;
                W = diag(w(rho_z).*dz);
                wTz = W*Tz;
                
                TzTz = (Tz.')*wTz;
                M = TT + TzTz;
                B = (T.')*rho + min_rho_z*transpose(sum(wTz));
                m_new = M\B;
                
                max_dm = max(abs(m_new-m)./max(abs(m_new)));
                iterations = iterations + 1;
                
                m = m_new;
            end
            
        end