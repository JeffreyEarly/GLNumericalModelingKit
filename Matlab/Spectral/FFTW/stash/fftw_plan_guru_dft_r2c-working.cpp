//
//  fftw_plan_guru_dft_r2c.cpp
//  
//
//  Created by Jeffrey Early on 12/10/24.
//

#include "mex.hpp"
#include "mexAdapter.hpp"
#include <fftw3.h>
#include <vector>

using matlab::mex::ArgumentList;
using matlab::data::TypedArray;
using matlab::data::ArrayFactory;
using matlab::data::ArrayDimensions;

class MexFunction : public matlab::mex::Function {
public:
    void operator()(ArgumentList outputs, ArgumentList inputs) {
        // Validate inputs
        if (inputs.size() != 2) {
            getEngine()->feval(u"error", 0,
                std::vector<matlab::data::Array>{ArrayFactory().createScalar("Two inputs required: matrix and dimension.")});
        }

        // Get the input matrix and transform dimension
        TypedArray<double> inputMatrix = inputs[0];
        int transformDim = static_cast<int>(inputs[1][0]);

        // Validate transform dimension
        ArrayDimensions dims = inputMatrix.getDimensions();
        if (transformDim < 1 || transformDim > static_cast<int>(dims.size())) {
            getEngine()->feval(u"error", 0,
                std::vector<matlab::data::Array>{ArrayFactory().createScalar("Invalid transform dimension.")});
        }

        // Convert MATLAB dimension (1-based) to C++ dimension (0-based)
        transformDim -= 1;

        // Determine the output dimensions (complex-valued along the transformed dimension)
        ArrayDimensions outputDims = dims;
        outputDims[transformDim] = dims[transformDim] / 2 + 1; // FFTW real-to-complex output size

        // Create output matrix
        ArrayFactory factory;
        TypedArray<std::complex<double>> outputMatrix = factory.createArray<std::complex<double>>(outputDims);

        // Prepare FFTW dimensions
        std::vector<fftw_iodim> fftDims(1); // Single dimension for the transform
        fftDims[0].n = dims[transformDim];
        fftDims[0].is = 1;
        fftDims[0].os = 1;

        std::vector<fftw_iodim> howmanyDims;
        for (size_t i = 0; i < dims.size(); ++i) {
            if (static_cast<int>(i) != transformDim) {
                fftw_iodim dim = { static_cast<int>(dims[i]), static_cast<int>(stride(dims, i)), static_cast<int>(stride(outputDims, i)) };
                howmanyDims.push_back(dim);
            }
        }

        double* fftwIn = &(*inputMatrix.begin());
        std::complex<double>* fftwOut = &(*outputMatrix.begin());
//        // Allocate FFTW buffers
//        std::vector<double> fftwIn(inputMatrix.getNumberOfElements());
//        std::vector<fftw_complex> fftwOut(outputMatrix.getNumberOfElements());
//
//        // Copy input data into FFTW input buffer
//        copyToFFTWBuffer(inputMatrix, fftwIn);

        // Create FFTW plan
        fftw_plan plan = fftw_plan_guru_dft_r2c(
            static_cast<int>(fftDims.size()), fftDims.data(),
            static_cast<int>(howmanyDims.size()), howmanyDims.data(),
            fftwIn, (fftw_complex *) fftwOut, FFTW_ESTIMATE);

        // Execute the FFT
        fftw_execute(plan);

        // Copy FFTW output buffer back to MATLAB array
//        copyFromFFTWBuffer(fftwOut, outputMatrix);

        // Clean up FFTW plan
        fftw_destroy_plan(plan);

        // Return the result
        outputs[0] = outputMatrix;
    }

private:
    size_t stride(const ArrayDimensions& dims, size_t dim) {
        size_t s = 1;
        for (size_t i = dim + 1; i < dims.size(); ++i) {
            s *= dims[i];
        }
        return s;
    }

    void copyToFFTWBuffer(const TypedArray<double>& inputMatrix, std::vector<double>& fftwIn) {
        auto it = inputMatrix.begin();
        for (size_t i = 0; i < fftwIn.size(); ++i, ++it) {
            fftwIn[i] = *it;
        }
    }

    void copyFromFFTWBuffer(const std::vector<fftw_complex>& fftwOut, TypedArray<std::complex<double>>& outputMatrix) {
        auto it = outputMatrix.begin();
        for (size_t i = 0; i < fftwOut.size(); ++i, ++it) {
            *it = std::complex<double>(fftwOut[i][0], fftwOut[i][1]);
        }
    }
};

