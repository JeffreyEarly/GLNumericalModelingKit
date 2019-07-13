
classdef BSplineUnitTests < matlab.unittest.TestCase

    methods (Test)
        function plusTest(testCase)
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.AbsoluteTolerance
            f = @(x) x;
            N = 11; % number of points
            K = 2; % order of spline
            
            % first let's do a uniform grid, lower order
            x = linspace(-1,1,N)';
            spline = InterpolatingSpline(x,f(x),K);
            
            spline = spline + 1;

            expSolution = f(x)+1;
            actSolution = spline(x);
            
            testCase.assertThat(actSolution, IsEqualTo(expSolution, 'Within', AbsoluteTolerance(2*eps))) 
        end
        
        function timesTest(testCase)
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.AbsoluteTolerance
            f = @(x) x;
            N = 11; % number of points
            K = 2; % order of spline
            
            % first let's do a uniform grid, lower order
            x = linspace(-1,1,N)';
            spline = InterpolatingSpline(x,f(x),K);
            
            spline = -2*spline;
            
            expSolution = -2*f(x);
            actSolution = spline(x);
            
            testCase.assertThat(actSolution, IsEqualTo(expSolution, 'Within', AbsoluteTolerance(2*eps))) 
        end
        
        function integrationTest(testCase)
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.AbsoluteTolerance
            f = @(x) x + ones(size(x));
            g = @(x) 0.5*x.^2 + x; % \int (x+1) dx = x^2/2 + x + const
            N = 11; % number of points
            K = 4; % order of spline
            
            % first let's do a uniform grid, lower order
            x = linspace(-1,1,N)';
            spline = InterpolatingSpline(x,f(x),K);
            
            intspline = cumsum(spline);
            
            expSolution = g(x)-g(x(1));
            actSolution = intspline(x);
            
            % There's a zero here, so we have to use absolute tolerance
            % only.
            testCase.assertThat(actSolution, IsEqualTo(expSolution, 'Within', AbsoluteTolerance(10*eps)))            
        end
        
        function differentiationTest(testCase)
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.RelativeTolerance
            f = @(x) -x.^3 + x.^2 - 2*x + 1;
            g = @(x) -3*x.^2 + 2*x - 2;
            N = 11; % number of points
            K = 4; % order of spline
            
            % first let's do a uniform grid, lower order
            x = linspace(-1,1,N)';
            spline = InterpolatingSpline(x,f(x),K);
            
            spline = diff(spline);
            
            expSolution = g(x);
            actSolution = spline(x);
            
            testCase.assertThat(actSolution, IsEqualTo(expSolution, 'Within', RelativeTolerance(100*eps))) 
        end
        
        function squareRootTest(testCase)
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.AbsoluteTolerance
            f = @(x) x.^2 - 2*x + 1;
            N = 11; % number of points
            K = 4; % order of spline
            
            % first let's do a uniform grid, lower order
            x = linspace(-1,1,N)';
            spline = InterpolatingSpline(x,f(x),K);
            
            spline = sqrt(spline);
            
            expSolution = sqrt(f(x));
            actSolution = spline(x);
            
            testCase.assertThat(actSolution, IsEqualTo(expSolution, 'Within', AbsoluteTolerance(2*eps))) 
        end
        
        
        function rootsTest(testCase)
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.AbsoluteTolerance
            f = @(x) mod(x,2)-0.5;
            N = 11; % number of points
            K = 2; % order of spline
            
            % first let's do a uniform grid, lower order
            x = linspace(0,10,N)';
            spline = InterpolatingSpline(x,f(x),K);
            
            r = roots(spline);
            
            expSolution = (0:9)' + 0.5;
            actSolution = r;
            
            testCase.assertThat(actSolution, IsEqualTo(expSolution, 'Within', AbsoluteTolerance(2*eps)))
        end
        
    end

end