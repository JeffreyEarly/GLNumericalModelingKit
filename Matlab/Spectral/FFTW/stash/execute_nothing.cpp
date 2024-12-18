//
//  execute_nothing.cpp
//  
//
//  Created by Jeffrey Early on 12/7/24.
//

#include "mex.hpp"
#include "mexAdapter.hpp"

class MexFunction : public matlab::mex::Function {
private:
    std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr;
    
public:
    void operator()(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs) override {
        namespace mat = matlab::data;
        mat::ArrayFactory factory;

        // Validate input arguments
        if (inputs.size() != 2 || outputs.size() != 1) {
            matlabPtr->feval(u"error", {factory.createScalar("Usage: output = execute_nothing(input)")});
            return;
        }
        mat::TypedArray<double> A = inputs[0];
        const mat::TypedArray<double> B = inputs[1];
        
        auto itA = A.begin();
        auto itB = B.begin();
        for (; itA != A.end(); ++itA, ++itB) {
            *itA = *itA + *itB;
        }
        
//        A[0] = static_cast<double>(1);
//        size_t numElements = A.getNumberOfElements();
//        for (size_t i = 0; i < numElements; ++i) {
//            A[i] = static_cast<double>(i + 1);
////            A[i] = A[i] + B[i]; // Example: Fill with 1, 2, 3, ...
//        }

        outputs[0] = A;
    }
};
