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
//        double* in;
//        double* out;
    };

    std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr;
    
public:
    void operator()(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs) override {
        namespace mat = matlab::data;
        mat::ArrayFactory factory;

        // Validate input arguments
        if (inputs.size() != 2 || outputs.size() != 1) {
            matlabPtr->feval(u"error", {factory.createScalar("Usage: output = execute_dct_plan_mex(plan, input)")});
            return;
        }

        // Retrieve the FFTW plan handle
        uint64_t handleValue = inputs[0][0];
        PlanHandle* handle = reinterpret_cast<PlanHandle*>(handleValue);

        // Retrieve the input data
//        mat::TypedArray<double> inputArray = std::move(inputs[1]);
        mat::TypedArray<double> inputArray = inputs[1];
//        size_t dataSize = inputArray.getNumberOfElements();

//        mxArray* mxInputArray = mat::getMATLABObject(inputArray);
//        mxArray* mxUnsharedArray = mxUnshareArray(mxInputArray);
        
        // Validate input size
//        size_t expectedSize = handle->dims.n;  // Example: adjust to match plan's dimension
//        if (dataSize != expectedSize) {
//            matlabPtr->feval(u"error", {factory.createScalar("Input size does not match plan dimensions.")});
//            return;
//        }

        // Copy input data to FFTW input buffer
//        std::copy(inputArray.begin(), inputArray.end(), handle->in);

        // Execute the FFTW plan
        // Getting the pointer causes the copy... probably for safety?
//        double* dataPtr = &(*inputArray.begin());
        auto itA = inputArray.begin();
//        mat::TypedArray<double> outputArray = factory.createArray<double>(inputArray.getDimensions());
//        fftw_execute_r2r(handle->plan, dataPtr, dataPtr);

        // Create MATLAB output array
//        mat::ArrayFactory factory;
//        mat::TypedArray<double> outputArray = factory.createArray<double>(mxUnsharedArray);
//        std::copy(handle->out, handle->out + dataSize, outputArray.begin());

        // Return the result
        outputs[0] = inputArray;
    }

//private:
//    struct PlanHandle {
//        fftw_plan plan;
//        double* in;
//        double* out;
//        fftw_iodim dims;
//    };
//
//    std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();
};

