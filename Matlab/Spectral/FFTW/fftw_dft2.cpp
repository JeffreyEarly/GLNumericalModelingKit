//  fftw_plan_guru_dft_r2c.cpp
//  Created by Jeffrey Early on 12/10/24.

// The most important thing about the whole exercise, is PREVENT MEMORY COPIES. It is so hard to *not* trigger a memory copy in Matlab, but if you do, then you will destroy your performance. The difference between getDataPtr<double>, and (*inputArray.begin()) --- is enormous!! The former can avoid a memory copy, the latter, cannot.
// Method 1: no memory copy, fast
//        TypedArray<double> inputArray = inputs[2];
//        const double* inPtr = getDataPtr<double>(inputArray);
// Method 2: no memory copy, requires known array sizes, fast
//        TypedArray<double> inputArray = std::move(inputs[2]);
//        double* inPtr = & inputArray[0][0].operator double &();
// Method 3: memory copy, slow, super simple
//        double* inPtr = &(*inputArray.begin());

#include "mex.hpp"
#include "mexAdapter.hpp"
#include <fftw3.h>
#include <vector>
#include <stdexcept>
#include <memory>
#include <numeric>

using matlab::mex::ArgumentList;
using matlab::data::TypedArray;
using matlab::data::ArrayFactory;
using matlab::data::ArrayDimensions;

class MexFunction : public matlab::mex::Function {
    
    matlab::data::ArrayFactory factory;
    
private:
    enum TransformType { TransformTypeDFTForward, TransformTypeDFTInverse, TransformTypeR2RForward, TransformTypeR2RInverse};
    
    struct DFTPlanHandle {
        fftw_plan planForward;
        fftw_plan planInverse;
        std::vector<size_t> realMatrixDims;
        std::vector<size_t> complexMatrixDims;
    };
    
    struct R2RPlanHandle {
        fftw_plan plan;
        std::vector<size_t> matrixDims;
    };

    std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();

    void flattenDimensionsAndTransformIndices(std::vector<size_t>& indims, std::vector<size_t>& transformDims) {
        for (auto& num : transformDims) {
            num -= 1; // Convert 1-based to 0-based indexing
        }

        std::vector<size_t> dims;
        for (size_t i = 0; i < indims.size(); ++i) {
            if (indims[i] != 1) {
                dims.push_back(indims[i]);
            } else {
                for (auto& tDim : transformDims) {
                    if (tDim == i) throw std::invalid_argument("Transform dimension has length 1.");
                    if (tDim > i) --tDim;
                }
            }
        }
        indims = dims;
    }

    void flattenDimensions(std::vector<size_t>& indims) {
        indims.erase(std::remove(indims.begin(), indims.end(), 1), indims.end());
    }

    ArrayDimensions outputDimensionsFromInputDimensions(const std::vector<size_t>& indims, const std::vector<size_t>& transformDims) {
        ArrayDimensions outputMatrixDims = indims;
        outputMatrixDims[transformDims.back() - 1] = indims[transformDims.back() - 1] / 2 + 1;
        return outputMatrixDims;
    }

    // The most important function: https://www.mathworks.com/matlabcentral/answers/1573713-c-mex-data-api-how-can-i-get-the-raw-pointer-of-an-array-input-without-copying-it
    template <typename T>
    const T* getDataPtr(matlab::data::Array arr) {
        const matlab::data::TypedArray<T> arr_t = arr;
        matlab::data::TypedIterator<const T> it(arr_t.begin());
        return it.operator->();
    }

    std::vector<size_t> typedArrayToSizeTVector(const matlab::data::TypedArray<double>& inputArray) {
        std::vector<size_t> dims;
        dims.reserve(inputArray.getNumberOfElements());
        for (const auto& val : inputArray) {
            dims.push_back(static_cast<size_t>(val));
        }
        return dims;
    }

    void fftwDimsSetup(const std::vector<size_t>& dims, const std::vector<size_t>& transformDims,
                       std::vector<fftw_iodim>& iodims, std::vector<fftw_iodim>& howmany_dims, TransformType type, double* scaleFactor = nullptr) {
        std::vector<size_t> outdims(dims);

        if (type == TransformTypeDFTForward) {
            outdims[transformDims.back()] = dims[transformDims.back()] / 2 + 1;
        } else if (type == TransformTypeDFTInverse)  {
            outdims[transformDims.back()] = 2 * (dims[transformDims.back()] - 1);
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

        for (size_t i = 0, t = 0; i < alldims.size(); ++i) {
            if (t < transformDims.size() && i == transformDims[t]) {
                iodims.push_back(alldims[i]);
                if (type == TransformTypeDFTForward && scaleFactor) *scaleFactor *= 1.0 / dims[i];
                ++t;
            } else {
                howmany_dims.push_back(alldims[i]);
            }
        }

        if (type == TransformTypeDFTInverse) {
            iodims.back().n = outdims[transformDims.back()];
        }
    }
    
    void createR2RPlan(ArgumentList outputs, ArgumentList inputs) {
        enum r2rArguments {dimArg = 1, transformDimArg = 2, nCoresArg = 3, kindArg = 4, plannerArg = 5};
        
        auto realMatrixDims = typedArrayToSizeTVector(inputs[dimArg]);
        auto transformDims = typedArrayToSizeTVector(inputs[transformDimArg]);
        std::vector<size_t> flatRealMatrixDims = realMatrixDims;
        flattenDimensionsAndTransformIndices(flatRealMatrixDims, transformDims);
        
        std::vector<fftw_iodim> iodims, howmany_dims;
        fftwDimsSetup(flatRealMatrixDims, transformDims, iodims, howmany_dims, TransformTypeR2RForward);
        
        int nCores = static_cast<int>(inputs[nCoresArg][0]);
        fftw_r2r_kind kind = (fftw_r2r_kind) static_cast<int>(inputs[kindArg][0]);
        unsigned planner = static_cast<unsigned>(inputs[plannerArg][0]);
        
        fftw_init_threads();
        fftw_plan_with_nthreads(nCores);
        
        int totalSize = std::accumulate(flatRealMatrixDims.begin(), flatRealMatrixDims.end(), 1, std::multiplies<>());
        double* in = fftw_alloc_real(totalSize);
        double* out = fftw_alloc_real(totalSize);
        fftw_plan plan = fftw_plan_guru_r2r(static_cast<int>(iodims.size()), iodims.data(), static_cast<int>(howmany_dims.size()), howmany_dims.data(), in, out, &kind, planner);
        fftw_free(in);
        fftw_free(out);
        
        // Wrap the plan and buffers in a struct
        R2RPlanHandle* handle = new R2RPlanHandle{plan, realMatrixDims};

        // Return the handle as an opaque pointer
        outputs[0] = factory.createScalar(reinterpret_cast<uint64_t>(handle));
    }
    
    void freeR2RPlan(ArgumentList outputs, ArgumentList inputs) {
        auto handle = reinterpret_cast<R2RPlanHandle*>(static_cast<uint64_t>(inputs[1][0]));
        if (handle) {
            fftw_destroy_plan(handle->plan);
            delete handle;
        }
    }
    
    void r2r(ArgumentList outputs, ArgumentList inputs) {
        auto handle = reinterpret_cast<R2RPlanHandle*>(static_cast<uint64_t>(inputs[1][0]));
        TypedArray<double> inputArray = std::move(inputs[2]);
        double* dataPtr = &(*inputArray.begin());
        fftw_execute_r2r(handle->plan, dataPtr, dataPtr);
        outputs[0] = inputArray;
    }
    
    void r2r_inout(ArgumentList outputs, ArgumentList inputs) {
        auto handle = reinterpret_cast<R2RPlanHandle*>(static_cast<uint64_t>(inputs[1][0]));
        TypedArray<double> inputArray = inputs[2];
        TypedArray<double> outputArray = std::move(inputs[3]); // Necessary! Prevents a memory copy
        auto inPtr = getDataPtr<double>(inputArray); // Necessary! Prevents a memory copy.
        double* outPtr = &(*outputArray.begin());
        fftw_execute_r2r(handle->plan, (double *) inPtr, outPtr);
        outputs[0] = outputArray;
    }
    
    void createDFTPlan(ArgumentList outputs, ArgumentList inputs) {
        auto realMatrixDims = typedArrayToSizeTVector(inputs[1]);
        auto transformDims = typedArrayToSizeTVector(inputs[2]);
        auto complexMatrixDims = outputDimensionsFromInputDimensions(realMatrixDims, transformDims);

        std::vector<size_t> flatRealMatrixDims = realMatrixDims;
        std::vector<size_t> flatComplexMatrixDims = complexMatrixDims;
        flattenDimensionsAndTransformIndices(flatRealMatrixDims, transformDims);
        flattenDimensions(flatComplexMatrixDims);

        std::vector<fftw_iodim> iodims, howmany_dims;
        double scaleFactor = 1.0;
        fftwDimsSetup(flatRealMatrixDims, transformDims, iodims, howmany_dims, TransformTypeDFTForward, &scaleFactor);
        
        int totalRealSize = std::accumulate(flatRealMatrixDims.begin(), flatRealMatrixDims.end(), 1, std::multiplies<>());
        int totalComplexSize = std::accumulate(flatComplexMatrixDims.begin(), flatComplexMatrixDims.end(), 1, std::multiplies<>());
        double* in = fftw_alloc_real(totalRealSize);
        fftw_complex* out = fftw_alloc_complex(totalComplexSize);
        
        int nCores = static_cast<int>(inputs[3][0]);
        fftw_init_threads();
        fftw_plan_with_nthreads(nCores);

        unsigned planner = static_cast<unsigned>(inputs[4][0]);
        fftw_plan planForward = fftw_plan_guru_dft_r2c(iodims.size(), iodims.data(), howmany_dims.size(), howmany_dims.data(), in, out, planner);
        fftwDimsSetup(flatComplexMatrixDims, transformDims, iodims, howmany_dims, TransformTypeDFTInverse);
        printDims(iodims);
        printDims(howmany_dims);
        fftw_plan planInverse = fftw_plan_guru_dft_c2r(iodims.size(), iodims.data(), howmany_dims.size(), howmany_dims.data(), out, in, planner);
        fftw_free(in);
        fftw_free(out);
        
        auto handle = new DFTPlanHandle{planForward, planInverse, realMatrixDims, complexMatrixDims};

        outputs[0] = factory.createScalar(reinterpret_cast<uint64_t>(handle));
        auto dimsArray = factory.createArray<double>({complexMatrixDims.size()});
        std::transform(complexMatrixDims.begin(), complexMatrixDims.end(), dimsArray.begin(), [](size_t dim) { return static_cast<double>(dim); });
        outputs[1] = dimsArray;
        outputs[2] = factory.createScalar(scaleFactor);
    }

    void freeDFTPlan(ArgumentList outputs, ArgumentList inputs) {
        auto handle = reinterpret_cast<DFTPlanHandle*>(static_cast<uint64_t>(inputs[1][0]));
        if (handle) {
            fftw_destroy_plan(handle->planForward);
            fftw_destroy_plan(handle->planInverse);
            delete handle;
        }
    }

    void r2c(ArgumentList outputs, ArgumentList inputs) {
        auto handle = reinterpret_cast<DFTPlanHandle*>(static_cast<uint64_t>(inputs[1][0]));
        TypedArray<double> inputArray = inputs[2];
        const double* inPtr = getDataPtr<double>(inputArray);
        auto outputArray = ArrayFactory().createArray<std::complex<double>>(handle->complexMatrixDims);
        std::complex<double> * outPtr = &(*outputArray.begin());
        fftw_execute_dft_r2c(handle->planForward, (double*) inPtr, (fftw_complex *) outPtr);
        outputs[0] = outputArray;
    }

    void c2r(ArgumentList outputs, ArgumentList inputs) {
        auto handle = reinterpret_cast<DFTPlanHandle*>(static_cast<uint64_t>(inputs[1][0]));
        TypedArray<std::complex<double>> inputArray = inputs[2];
        auto outputArray = ArrayFactory().createArray<double>(handle->realMatrixDims);
        fftw_execute_dft_c2r(handle->planInverse, reinterpret_cast<fftw_complex*>(inputArray.begin().operator->()), outputArray.begin().operator->());
        outputs[0] = outputArray;
    }
    
    void r2c_inout(ArgumentList outputs, ArgumentList inputs) {
        auto handle = reinterpret_cast<DFTPlanHandle*>(static_cast<uint64_t>(inputs[1][0]));
        TypedArray<double> inputArray = inputs[2];
        const double* inPtr = getDataPtr<double>(inputArray);
        TypedArray<std::complex<double>> outputArray = std::move(inputs[3]); // Necessary! Prevents a memory copy
        const std::complex<double> * outPtr = getDataPtr<std::complex<double>>(outputArray); // Also necessary! Prevents a memory copy.
        fftw_execute_dft_r2c(handle->planForward, (double *) inPtr, (fftw_complex *) outPtr);
        outputs[0] = outputArray;
    }

public:
    void operator()(ArgumentList outputs, ArgumentList inputs) override {
        if (inputs[0].getType() == matlab::data::ArrayType::CHAR) {
            matlab::data::CharArray commandArray = inputs[0];
            auto command = commandArray.toAscii();
            if (command == "create") createDFTPlan(outputs, inputs);
            else if (command == "free") freeDFTPlan(outputs, inputs);
            else if (command == "r2c") r2c(outputs, inputs);
            else if (command == "c2r") c2r(outputs, inputs);
            else if (command == "r2c_inout") r2c_inout(outputs, inputs);
            else if (command == "createR2RPlan") createR2RPlan(outputs, inputs);
            else if (command == "r2r") r2r(outputs, inputs);
            else if (command == "r2r_inout") r2r_inout(outputs, inputs);
            else matlabPtr->feval(u"error", 0, {ArrayFactory().createScalar("Unknown command.")});
        } else {
            matlabPtr->feval(u"error", 0, {ArrayFactory().createScalar("Invalid input.")});
        }
    }
    
    void printVec(const std::vector<size_t>& alldims) {
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
    
    void printDims(const std::vector<fftw_iodim>& alldims) {
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
};
