//
//  mexPlan.c
//  
//
//  Created by Jeffrey Early on 12/5/24.
//

#include <complex>
#include "fftw3.h"
#include "mex.hpp"
#include "mexAdapter.hpp"

class MexFunction : public matlab::mex::Function {
public:
    void operator()(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs) override {
        namespace mat = matlab::data;
        mat::ArrayFactory factory;

        // Check number of inputs/outputs
        if (inputs.size() != 1 || outputs.size() != 1) {
            matlabPtr->feval(u"error", {factory.createScalar("One input and one output required.")});
            return;
        }

        // Extract input (assumed to be a 1D double array)
        mat::TypedArray<double> inputArray = inputs[0];
        size_t n = inputArray.getNumberOfElements();

        // Allocate FFTW arrays
        fftw_complex *in = fftw_alloc_complex(n);
        fftw_complex *out = fftw_alloc_complex(n);

        // Copy input data into FFTW input array
        for (size_t i = 0; i < n; ++i) {
            in[i][0] = inputArray[i]; // Real part
            in[i][1] = 0.0;          // Imaginary part
        }

        // Plan and execute FFT
        fftw_plan plan = fftw_plan_dft_1d(n, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
        fftw_execute(plan);

        // Create MATLAB output array
        mat::TypedArray<std::complex<double>> outputArray = factory.createArray<std::complex<double>>({n});
        for (size_t i = 0; i < n; ++i) {
            outputArray[i] = std::complex<double>(out[i][0], out[i][1]);
        }

        // Set the output
        outputs[0] = outputArray;

        // Cleanup
        fftw_destroy_plan(plan);
        fftw_free(in);
        fftw_free(out);
    }

private:
    std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();
};
