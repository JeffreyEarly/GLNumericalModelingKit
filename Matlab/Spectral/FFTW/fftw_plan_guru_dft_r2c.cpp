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
private:
    struct PlanHandle {
        fftw_plan plan;
        std::vector<size_t> inputMatrixDims;
        std::vector<size_t> outputMatrixDims;
    };

//    std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr;
    
    ArrayDimensions createFFTWIODim(const std::vector<size_t>& indims, std::vector<size_t>& transformDims, std::vector<fftw_iodim>& iodims, std::vector<fftw_iodim>& howmany_dims) {
        // Validate transformDim
//        if (transformDims < 1 || transformDims > static_cast<int>(indims.size())) {
//            throw std::invalid_argument("Invalid transform dimension");
//        }
        
        // Convert 1-based indexing to 0-based indexing
        for (size_t &num : transformDims) {
            num -= 1;
        }
        
        // Determine the output dimensions (complex-valued along the transformed dimension)
        ArrayDimensions outputMatrixDims = indims;
        outputMatrixDims[transformDims.back()] = indims[transformDims.back()] / 2 + 1;
        
        // imagine indims={10, 1, 5, 1, 2} and transformDims={0, 2}, this should become
        // dims={10,5,2} and transformDims={0,1}. According to what I have read online, FFTW
        // is fine with singleton dimensions, but other FFT libraries are not.
        std::vector<size_t> dims;
        for (size_t i = 0; i < indims.size(); i++)  {
            if (indims[i] == 1) {
                for (size_t &tDim : transformDims) {
                    if (tDim > i) {
                        tDim -= 1;
                    } else if ( tDim == i) {
                        throw std::invalid_argument("Transform dimension has length 1.");
                    }
                }
            } else {
                dims.push_back(indims[i]);
            }
        }
        
        std::vector<size_t> outdims;
        for (size_t i = 0; i < dims.size(); i++) {
            if (transformDims.back() == i) {
                outdims.push_back(dims[i]/2 + 1);
            } else {
                outdims.push_back(dims[i]);
            }
        }
        
        std::vector<fftw_iodim> alldims;
        for (size_t i = 0; i < dims.size(); i++) {
            fftw_iodim newdim;
            newdim.n = dims[i];
            if (alldims.size() == 0) {
                newdim.is = 1;
                newdim.os = 1;
            } else {
                size_t n = alldims.size();
                newdim.is = alldims[n-1].is*dims[n-1];
                newdim.os = alldims[n-1].os*outdims[n-1];
            }
            alldims.push_back(newdim);
        }
        matlab::data::ArrayFactory factory;
        std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();
        for (size_t i = 0; i < alldims.size(); i++) {
            matlabPtr->feval(u"fprintf", 0,
                     std::vector<matlab::data::Array>(
                      {factory.createScalar("alldims (%d, %d, %d)\n"),
                          factory.createScalar(alldims[i].n),
                          factory.createScalar(alldims[i].is),
                          factory.createScalar(alldims[i].os)}));
        }
        
        iodims.clear();
        howmany_dims.clear();
        size_t iTransformDim = 0;
        for (size_t i = 0; i < alldims.size(); ++i) {
            if (iTransformDim < transformDims.size() && i == transformDims[iTransformDim]) {
                iodims.push_back(alldims[i]);
                iTransformDim += 1;
            } else {
                howmany_dims.push_back(alldims[i]);
                
            }
        }
        

//        matlabPtr->feval(u"fprintf", 0, std::vector<matlab::data::Array>( {factory.createScalar("iodims (%i, %i, %i)\n")}));
        for (size_t i = 0; i < iodims.size(); i++) {
            matlabPtr->feval(u"fprintf", 0,
                     std::vector<matlab::data::Array>(
                      {factory.createScalar("iodims (%d, %d, %d)\n"),
                          factory.createScalar(iodims[i].n),
                          factory.createScalar(iodims[i].is),
                          factory.createScalar(iodims[i].os)}));
        }
        
        for (size_t i = 0; i < howmany_dims.size(); i++) {
            matlabPtr->feval(u"fprintf", 0,
                     std::vector<matlab::data::Array>(
                      {factory.createScalar("howmany_dims (%d, %d, %d)\n"),
                          factory.createScalar(howmany_dims[i].n),
                          factory.createScalar(howmany_dims[i].is),
                          factory.createScalar(howmany_dims[i].os)}));
        }
        
        return outputMatrixDims;
    }
    
    std::vector<size_t> typedArrayToSizeTVector(const matlab::data::TypedArray<double>& inputArray) {
        std::vector<size_t> dims(inputArray.getNumberOfElements());
        size_t index = 0;
        for (const auto& val : inputArray) {
            dims[index++] = static_cast<size_t>(val);
        }
        return dims;
    }
public:
    void operator()(ArgumentList outputs, ArgumentList inputs) {
        // Validate inputs
//        if (inputs.size() != 2) {
//            getEngine()->feval(u"error", 0,
//                std::vector<matlab::data::Array>{ArrayFactory().createScalar("Two inputs required: matrix and dimension.")});
//        }

        std::vector<size_t> inputDims = typedArrayToSizeTVector(inputs[0]);
        std::vector<size_t> transformDims = typedArrayToSizeTVector(inputs[1]);
        std::vector<fftw_iodim> iodims;
        std::vector<fftw_iodim> howmany_dims;
        ArrayDimensions outputDims = createFFTWIODim(inputDims, transformDims, iodims, howmany_dims);
        
        // Create output matrix
        ArrayFactory factory;
        TypedArray<std::complex<double>> outputMatrix = factory.createArray<std::complex<double>>(outputDims);

        int totalSize = 1;
        for (const auto& dim : inputDims) {
            totalSize *= dim;
        }
        double* in = fftw_alloc_real(totalSize);
        
        totalSize = 1;
        for (const auto& dim : outputDims) {
            totalSize *= dim;
        }
        fftw_complex* out = fftw_alloc_complex(totalSize);

        // Create FFTW plan
        fftw_plan plan = fftw_plan_guru_dft_r2c(
            static_cast<int>(iodims.size()), iodims.data(),
            static_cast<int>(howmany_dims.size()), howmany_dims.data(),
            in, out, FFTW_ESTIMATE);
        fftw_free(in);
        fftw_free(out);
        
        // Wrap the plan and buffers in a struct
        PlanHandle* handle = new PlanHandle{plan, inputDims, outputDims};

        // Return the handle as an opaque pointer
        outputs[0] = factory.createScalar(reinterpret_cast<uint64_t>(handle));
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

