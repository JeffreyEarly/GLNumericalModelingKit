classdef StudentTDistribution < Distribution
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        sigma
        nu
    end
    
    methods
        function self = StudentTDistribution(sigma,nu)
            self.sigma = sigma;
            self.nu = nu;
            self.pdf = @(z) gamma((nu+1)/2)./(sqrt(pi*nu)*sigma*gamma(nu/2)*(1+(z.*z)/(nu*sigma*sigma)).^((nu+1)/2));
            self.cdf = @(z) StudentTDistribution.tcdf(z/sigma,nu);
            self.w = @(z)((nu/(nu+1))*sigma^2*(1+z.^2/(nu*sigma^2)));
            self.variance = sigma*sigma*nu/(nu-2);
            
            a = gamma((nu+1)/2)./(sqrt(pi*nu)*sigma*gamma(nu/2));
            c = nu*sigma*sigma;
            m = (nu+1)/2;
            % so that
            % pdf = gamma(m)./(sqrt(c)*gamma(m-1/2)*(1+(z.*z)/c).^(m));
            % log(pdf) = - m*log(1+(z.*z)/c) + log(gamma(m))-log(gamma(m-1/2)) - log(sqrt(c))  
            self.dPDFoverZ = @(z) -2*(a*m/c)*(1+z.*z/c).^(-m-1);
            self.logPDF = @(z) -m*log(1+z.*z/c) + length(z)*(log(gamma(m))-log(sqrt(pi*nu)*sigma*gamma(nu/2)));
            % ln \Gamma(z) = \ln \Gamma(z+1) - \ln z
            % ln \Gamma(z) = \ln \Gamma(z+1) - \ln z
        end
        
    end
    
    methods (Static)
       function p = tcdf(x,n)
            % TCDF returns student cumulative distribtion function
            %
            % cdf = tcdf(x,DF);
            %
            % Computes the CDF of the students distribution
            %    with DF degrees of freedom
            % x,DF must be matrices of same size, or any one can be a scalar.
            %
            % see also: NORMCDF, NORMPDF, NORMINV
            
            % Reference(s):
            
            %	$Revision: 1.1 $
            %	$Id: tcdf.m,v 1.1 2003/09/12 12:14:45 schloegl Exp $
            %	Copyright (c) 2000-2003 by  Alois Schloegl <a.schloegl@ieee.org>
            
            %    This program is free software; you can redistribute it and/or modify
            %    it under the terms of the GNU General Public License as published by
            %    the Free Software Foundation; either version 2 of the License, or
            %    (at your option) any later version.
            %
            %    This program is distributed in the hope that it will be useful,
            %    but WITHOUT ANY WARRANTY; without even the implied warranty of
            %    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
            %    GNU General Public License for more details.
            %
            %    You should have received a copy of the GNU General Public License
            %    along with this program; if not, write to the Free Software
            %    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
            
            
            % check size of arguments
            n = x+n-x;	  % if this line causes an error, size of input arguments do not fit.
            z = n ./ (n + x.^2);
            
            % allocate memory
            p = z;
            p(x==Inf) = 1;
            
            % workaround for invalid arguments in BETAINC
            tmp   = isnan(z) | ~(n>0);
            p(tmp)= NaN;
            ix    = (~tmp);
            p(ix) = betainc (z(ix), n(ix)/2, 1/2) / 2;
            
            ix    = find(x>0);
            p(ix) = 1 - p(ix);
            
            % shape output
            p = reshape(p,size(z));
        end 
    end
end

