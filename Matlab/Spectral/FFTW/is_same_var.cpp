#include "mex.hpp"
#include "mexAdapter.hpp"
#include <vector>
#include <stdexcept>
#include <memory>
#include <numeric>

using matlab::mex::ArgumentList;
using matlab::data::TypedArray;
using matlab::data::ArrayFactory;
using matlab::data::ArrayDimensions;

class MexFunction : public matlab::mex::Function {
    
    template <typename T>
    const T* getDataPtr(matlab::data::Array arr) {
        const matlab::data::TypedArray<T> arr_t = arr;
        matlab::data::TypedIterator<const T> it(arr_t.begin());
        return it.operator->();
    }
    
   

    
public:
    void operator()(ArgumentList outputs, ArgumentList inputs) override {
        
        std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();
        
        TypedArray<double> aArray = inputs[0];
        TypedArray<double> bArray = inputs[1];

//        auto aPtr = &(*aArray.begin());
//        auto bPtr = &(*bArray.begin());
        
        auto aPtr = getDataPtr<double>(aArray);
        auto bPtr = getDataPtr<double>(bArray);
        
        matlab::data::ArrayFactory factory;
//        matlabPtr->feval(u"fprintf", 0, std::vector<matlab::data::Array>({factory.createScalar("address 0x%x\n"),
//                                                                             factory.createScalar(aPtr)}));
//        
//        matlabPtr->feval(u"fprintf", 0,
//                         std::vector<matlab::data::Array>({factory.createScalar(aPtr)});
////            matlabPtr->feval(u"fprintf", 0,
//                     std::vector<matlab::data::Array>({factory.createScalar(bPtr}));

        // Check if the pointers are the same
        bool isSameVariable = (aPtr == bPtr);

        // Return the result to MATLAB as a logical
        matlab::data::ArrayFactory factory;
        outputs[0] = factory.createScalar(isSameVariable);
    }
};
