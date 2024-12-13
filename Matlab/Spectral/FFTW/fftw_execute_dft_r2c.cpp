//
//  execute_dct_plan_mex.cpp
//  
//
//  Created by Jeffrey Early on 12/7/24.
//

#include "mex.hpp"
#include "mexAdapter.hpp"
#include <fftw3.h>

using namespace matlab::mex;
using namespace matlab::data;

//extern "C" bool mxUnshareArray(mxArray *array_ptr, bool noDeepCopy);

//extern "C" {
//    #include "matrix.h" // For mxUnshareArray
//}

class MexFunction : public matlab::mex::Function {
private:
    struct PlanHandle {
        fftw_plan plan;
        std::vector<size_t> inputMatrixDims;
        std::vector<size_t> outputMatrixDims;
    };
    std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr;
    
    template <typename T>
    const T* getDataPtr(matlab::data::Array arr) {
      const matlab::data::TypedArray<T> arr_t = arr;
      matlab::data::TypedIterator<const T> it(arr_t.begin());
      return it.operator->();
    }
    
public:
    void operator()(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs) override {
        namespace mat = matlab::data;
        mat::ArrayFactory factory;

        // Validate input arguments
        if (inputs.size() != 2 || outputs.size() != 1) {
            matlabPtr->feval(u"error", {factory.createScalar("Usage: output = fftw_execute_dft_r2c(plan, input)")});
            return;
        }
        
        // Retrieve the FFTW plan handle
        uint64_t handleValue = inputs[0][0];
        PlanHandle* handle = reinterpret_cast<PlanHandle*>(handleValue);
        
        mat::TypedArray<double> inputArray = inputs[1];
        TypedArray<std::complex<double>> outputArray = factory.createArray<std::complex<double>>(handle->outputMatrixDims);

        auto inPtr = getDataPtr<double>(inputArray);
        std::complex<double> * outPtr = &(*outputArray.begin());
        
        fftw_execute_dft_r2c(handle->plan, (double *) inPtr, (fftw_complex*) outPtr);
        outputs[0] = outputArray;
                
//        mat::TypedArray<double> inputArray = inputs[1];
//        mat::TypedArray<double> outputArray = std::move(inputs[2]);
//
////        double* inPtr = &(*inputArray.begin());
//        auto inPtr = getDataPtr<double>(inputArray);
//        double* outPtr = &(*outputArray.begin());
////        auto outPtr = getDataPtr<double>(outputArray);
//
//        fftw_execute_r2r(handle->plan, (double *) inPtr, outPtr);
//        outputs[0] = outputArray;
    }
    

};

