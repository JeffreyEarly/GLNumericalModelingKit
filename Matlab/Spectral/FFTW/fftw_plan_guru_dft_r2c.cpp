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
        fftw_plan planForward;
        fftw_plan planInverse;
        std::vector<size_t> realMatrixDims;
        std::vector<size_t> complexMatrixDims;
    };

    std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr;
    
    // Removes the singleton dimensions from indims, and adjust the transformDims as appropriate.
    void flattenDimensionsAndTransformIndices(std::vector<size_t>& indims, std::vector<size_t>& transformDims)
    {
        // Convert 1-based indexing to 0-based indexing
        for (size_t &num : transformDims) {
            num -= 1;
        }
        
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
        indims = dims;
    }
    
    void flattenDimensions(std::vector<size_t>& indims)
    {
        std::vector<size_t> dims;
        for (size_t i = 0; i < indims.size(); i++)  {
            if (indims[i] > 1) {
                dims.push_back(indims[i]);
            }
        }
        indims = dims;
    }
    
    // Returns the appropriate output dimensions for a complex matrix in matlab *with* singletons intact.
    // This assume indims is the un-altered (not flattened) dimensions of the input matrix.
    ArrayDimensions outputDimensionsFromInputDimensions(const std::vector<size_t>& indims, const std::vector<size_t>& transformDims)
    {
        ArrayDimensions outputMatrixDims = indims;
        outputMatrixDims[transformDims.back()-1] = indims[transformDims.back()-1] / 2 + 1;
        return outputMatrixDims;
    }
    
    void fftwForwardDims(const std::vector<size_t>& dims, const std::vector<size_t>& transformDims, std::vector<fftw_iodim>& iodims, std::vector<fftw_iodim>& howmany_dims, double *scaleFactor)
    {
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

        iodims.clear();
        howmany_dims.clear();
        size_t iTransformDim = 0;
        for (size_t i = 0; i < alldims.size(); ++i) {
            if (iTransformDim < transformDims.size() && i == transformDims[iTransformDim]) {
                iodims.push_back(alldims[i]);
                iTransformDim += 1;
                *scaleFactor *= 1.0 / ((double) alldims[i].n);
            } else {
                howmany_dims.push_back(alldims[i]);
                
            }
        }
    }
 
    void fftwInverseDims(const std::vector<size_t>& dims, const std::vector<size_t>& transformDims, std::vector<fftw_iodim>& iodims, std::vector<fftw_iodim>& howmany_dims)
    {
        std::vector<size_t> outdims;
        for (size_t i = 0; i < dims.size(); i++) {
            if (transformDims.back() == i) {
                outdims.push_back(2*(dims[i] - 1));
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
        
        iodims.back().n = 2*(iodims.back().n-1);
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
        matlab::data::ArrayFactory factory;
        std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();
        
        std::vector<size_t> realMatrixDims = typedArrayToSizeTVector(inputs[0]);
        std::vector<size_t> transformDims = typedArrayToSizeTVector(inputs[1]); // *index* of the dim to transform
        ArrayDimensions complexMatrixDims = outputDimensionsFromInputDimensions(realMatrixDims, transformDims);
        
        std::vector<size_t> flatRealMatrixDims = realMatrixDims;
        flattenDimensionsAndTransformIndices(flatRealMatrixDims, transformDims);
        
        std::vector<size_t> flatComplexMatrixDims = complexMatrixDims;
        flattenDimensions(flatComplexMatrixDims);
        
        std::vector<fftw_iodim> iodims;
        std::vector<fftw_iodim> howmany_dims;
        double scaleFactor = 1.0;
        fftwForwardDims(flatRealMatrixDims, transformDims, iodims, howmany_dims, &scaleFactor);
        printDims(iodims);
        printDims(howmany_dims);
        printVec(flatRealMatrixDims);
        printVec(flatComplexMatrixDims);
        
        //createFFTWIODim(realMatrixDims, transformDims, iodims, howmany_dims, &scaleFactor);
        
        int totalSize = 1;
        matlabPtr->feval(u"fprintf", 0, std::vector<matlab::data::Array>({factory.createScalar("totalSize: %d\n"), factory.createScalar(totalSize)}));
        
        for (const auto& dim : flatRealMatrixDims) {
            totalSize *= dim;
        }
        matlabPtr->feval(u"fprintf", 0, std::vector<matlab::data::Array>({factory.createScalar("totalSize: %d\n"), factory.createScalar(totalSize)}));
        double* in = fftw_alloc_real(totalSize);
        
        totalSize = 1;
        for (const auto& dim : flatComplexMatrixDims) {
            totalSize *= dim;
        }
        matlabPtr->feval(u"fprintf", 0, std::vector<matlab::data::Array>({factory.createScalar("totalSize: %d\n"), factory.createScalar(totalSize)}));
        fftw_complex* out = fftw_alloc_complex(totalSize);
        


        // Create FFTW plan
        int nCores = static_cast<int>(inputs[2][0]);
        fftw_init_threads();
        fftw_plan_with_nthreads(nCores);
        
        unsigned planner = static_cast<unsigned>(inputs[3][0]);
        fftw_plan planForward = fftw_plan_guru_dft_r2c(
            static_cast<int>(iodims.size()), iodims.data(),
            static_cast<int>(howmany_dims.size()), howmany_dims.data(),
            in, out, planner);
        
        std::vector<fftw_iodim> iodims_inverse;
        std::vector<fftw_iodim> howmany_dims_inverse;
        fftwInverseDims(complexMatrixDims, transformDims, iodims_inverse, howmany_dims_inverse);
        printDims(iodims_inverse);
        printDims(howmany_dims_inverse);
        
        fftw_plan planInverse = fftw_plan_guru_dft_c2r(
            static_cast<int>(iodims_inverse.size()), iodims_inverse.data(),
            static_cast<int>(howmany_dims_inverse.size()), howmany_dims_inverse.data(),
            out, in, planner);
        
        fftw_free(in);
        fftw_free(out);
        
        
        
        // Wrap the plan and buffers in a struct
        PlanHandle* handle = new PlanHandle{planForward, planInverse, realMatrixDims, complexMatrixDims};

        // Return the handle as an opaque pointer
//        ArrayFactory factory;
        outputs[0] = factory.createScalar(reinterpret_cast<uint64_t>(handle));
        
        // Return the size of the complex output array
        TypedArray<double> dimsArray = factory.createArray<double>({complexMatrixDims.size()});
        size_t index = 0;
        for (const auto& dim : complexMatrixDims) {
            dimsArray[index++] = static_cast<double>(dim); // MATLAB uses doubles for numeric values
        }
        outputs[1] = dimsArray;
        
        // Return the scale factor for the forward transform
        outputs[2] = factory.createScalar(scaleFactor);
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
    
//    void print(string message)
//    {
//        
//    }
    
    void printDims(const std::vector<fftw_iodim>& alldims) {
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
    }
    
    void printVec(const std::vector<size_t>& alldims) {
        matlab::data::ArrayFactory factory;
        std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();
        matlabPtr->feval(u"fprintf", 0,
                 std::vector<matlab::data::Array>({factory.createScalar("[")}));
        for (size_t i = 0; i < alldims.size(); i++) {
            matlabPtr->feval(u"fprintf", 0,
                     std::vector<matlab::data::Array>(
                      {factory.createScalar(" %d"),
                          factory.createScalar(alldims[i])}));
        }
        matlabPtr->feval(u"fprintf", 0,
                 std::vector<matlab::data::Array>({factory.createScalar("]\n")}));
    }
};

// Useful debugging code
//matlab::data::ArrayFactory factory;
//std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();
//for (size_t i = 0; i < alldims.size(); i++) {
//    matlabPtr->feval(u"fprintf", 0,
//             std::vector<matlab::data::Array>(
//              {factory.createScalar("alldims (%d, %d, %d)\n"),
//                  factory.createScalar(alldims[i].n),
//                  factory.createScalar(alldims[i].is),
//                  factory.createScalar(alldims[i].os)}));
//}
