//
//  removeNegativeNumbers.cpp
//  
//
//  Created by Jeffrey Early on 12/8/24.
//

#include "mex.hpp"
#include "mexAdapter.hpp"

using namespace matlab::data;
using matlab::mex::ArgumentList;

class MexFunction : public matlab::mex::Function {
public:
    void operator()(ArgumentList outputs, ArgumentList inputs) {
        checkArguments(inputs);
        TypedArray<double> largeArray = std::move(inputs[0]);
//        auto elem = largeArray.begin();
        for (auto& elem : largeArray) {
////            elem = elem + 1;
            if (elem < 0) {
                elem = 0;
            }
        }
        outputs[0] = largeArray;
    }

    void checkArguments(ArgumentList inputs) {
        std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();
        ArrayFactory factory;

        if (inputs[0].getType() != ArrayType::DOUBLE ||
            inputs[0].getType() == ArrayType::COMPLEX_DOUBLE) {
            matlabPtr->feval(u"error", 0,
                std::vector<Array>({ factory.createScalar("Incorrect input") }));
        }
    }
};
