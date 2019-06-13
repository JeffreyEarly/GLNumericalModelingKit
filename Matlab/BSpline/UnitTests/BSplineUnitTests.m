classdef BSplineUnitTests < matlab.unittest.TestCase

    methods (Test)
        function testAddition(testCase)
            
            
            actSolution = quadraticSolver(1,-3,2);
            expSolution = [2,1];
            testCase.verifyEqual(actSolution,expSolution)
        end
    end

end