classdef TestableFFTWBackend < FFTWBackend
    methods (Static)
        function context = context()
            context = FFTWBackend.defaultContext;
        end

        function capabilities = inspect(context)
            capabilities = FFTWBackend.capabilitiesWithContext(context);
        end

        function capabilities = buildUsing(context)
            capabilities = FFTWBackend.buildWithContext(context);
        end
    end
end
