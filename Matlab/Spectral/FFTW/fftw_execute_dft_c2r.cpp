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

class MexFunction : public matlab::mex::Function {
private:
    struct PlanHandle {
        fftw_plan planForward;
        fftw_plan planInverse;
        std::vector<size_t> realMatrixDims;
        std::vector<size_t> complexMatrixDims;
    };
    std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();
    
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
            matlabPtr->feval(u"error", 0, std::vector<matlab::data::Array>({factory.createScalar("Usage: realMatrix = fftw_execute_dft_r2c(plan, complexMatrix)")}));
        }
        
        // Retrieve the FFTW plan handle
        uint64_t handleValue = inputs[0][0];
        PlanHandle* handle = reinterpret_cast<PlanHandle*>(handleValue);
        if (!handle) {
            matlabPtr->feval(u"error", 0, std::vector<matlab::data::Array>({factory.createScalar("Invalid plan provided")}));
        }
        
        if (inputs[1].getType() != matlab::data::ArrayType::COMPLEX_DOUBLE) {
            matlabPtr->feval(u"error", 0, std::vector<matlab::data::Array>({ factory.createScalar("Input matrix must be a complex double") }));
        }
        TypedArray<std::complex<double>> inputArray = inputs[1];
        
        if (inputArray.getDimensions() != handle->complexMatrixDims) {
            matlabPtr->feval(u"error", 0, std::vector<matlab::data::Array>({factory.createScalar("Input matrix dimensions do not match expected dimensions")}));
        }
        
        TypedArray<double> outputArray = factory.createArray<double>(handle->realMatrixDims);

        auto inPtr = getDataPtr<std::complex<double>>(inputArray);
        double * outPtr = &(*outputArray.begin());
        
        fftw_execute_dft_c2r(handle->planInverse, (fftw_complex*) inPtr, (double *) outPtr);
        outputs[0] = outputArray;
    }
};

