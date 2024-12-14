//  fftw_plan_guru_dft_r2c.cpp
//  Created by Jeffrey Early on 12/10/24.

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
private:
    struct PlanHandle {
        fftw_plan planForward;
        fftw_plan planInverse;
        std::vector<size_t> realMatrixDims;
        std::vector<size_t> complexMatrixDims;
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
                       std::vector<fftw_iodim>& iodims, std::vector<fftw_iodim>& howmany_dims, bool forward, double* scaleFactor = nullptr) {
        size_t scale = 1;
        std::vector<size_t> outdims(dims);

        if (forward) {
            outdims[transformDims.back()] = dims[transformDims.back()] / 2 + 1;
        } else {
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
                if (forward && scaleFactor) *scaleFactor *= 1.0 / dims[i];
                ++t;
            } else {
                howmany_dims.push_back(alldims[i]);
            }
        }

        if (!forward) iodims.back().n = outdims[transformDims.back()];
    }

    void createPlan(ArgumentList outputs, ArgumentList inputs) {
        auto realMatrixDims = typedArrayToSizeTVector(inputs[1]);
        auto transformDims = typedArrayToSizeTVector(inputs[2]);
        auto complexMatrixDims = outputDimensionsFromInputDimensions(realMatrixDims, transformDims);

        std::vector<size_t> flatRealMatrixDims = realMatrixDims;
        std::vector<size_t> flatComplexMatrixDims = complexMatrixDims;
        flattenDimensionsAndTransformIndices(flatRealMatrixDims, transformDims);
        flattenDimensions(flatComplexMatrixDims);

        std::vector<fftw_iodim> iodims, howmany_dims;
        double scaleFactor = 1.0;
        fftwDimsSetup(flatRealMatrixDims, transformDims, iodims, howmany_dims, true, &scaleFactor);
        printDims(iodims);
        printDims(howmany_dims);
        printVec(flatRealMatrixDims);
        printVec(flatComplexMatrixDims);
        
        int totalRealSize = std::accumulate(flatRealMatrixDims.begin(), flatRealMatrixDims.end(), 1, std::multiplies<>());
        int totalComplexSize = std::accumulate(flatComplexMatrixDims.begin(), flatComplexMatrixDims.end(), 1, std::multiplies<>());
        double* in = fftw_alloc_real(totalRealSize);
        fftw_complex* out = fftw_alloc_complex(totalComplexSize);
        
        int nCores = static_cast<int>(inputs[3][0]);
        fftw_init_threads();
        fftw_plan_with_nthreads(nCores);

        unsigned planner = static_cast<unsigned>(inputs[4][0]);
        fftw_plan planForward = fftw_plan_guru_dft_r2c(iodims.size(), iodims.data(), howmany_dims.size(), howmany_dims.data(), in, out, planner);
        fftwDimsSetup(flatComplexMatrixDims, transformDims, iodims, howmany_dims, false);
        printDims(iodims);
        printDims(howmany_dims);
        fftw_plan planInverse = fftw_plan_guru_dft_c2r(iodims.size(), iodims.data(), howmany_dims.size(), howmany_dims.data(), out, in, planner);
        fftw_free(in);
        fftw_free(out);
        
        auto handle = new PlanHandle{planForward, planInverse, realMatrixDims, complexMatrixDims};
        ArrayFactory factory;

        outputs[0] = factory.createScalar(reinterpret_cast<uint64_t>(handle));
        auto dimsArray = factory.createArray<double>({complexMatrixDims.size()});
        std::transform(complexMatrixDims.begin(), complexMatrixDims.end(), dimsArray.begin(), [](size_t dim) { return static_cast<double>(dim); });
        outputs[1] = dimsArray;
        outputs[2] = factory.createScalar(scaleFactor);
    }

    void freePlan(ArgumentList outputs, ArgumentList inputs) {
        auto handle = reinterpret_cast<PlanHandle*>(static_cast<uint64_t>(inputs[0][0]));
        if (handle) {
            fftw_destroy_plan(handle->planForward);
            fftw_destroy_plan(handle->planInverse);
            delete handle;
        }
    }

    void r2c(ArgumentList outputs, ArgumentList inputs) {
        auto handle = reinterpret_cast<PlanHandle*>(static_cast<uint64_t>(inputs[1][0]));
        TypedArray<double> inputArray = inputs[2];
        auto outputArray = ArrayFactory().createArray<std::complex<double>>(handle->complexMatrixDims);
        fftw_execute_dft_r2c(handle->planForward, const_cast<double*>(inputArray.begin().operator->()), reinterpret_cast<fftw_complex*>(outputArray.begin().operator->()));
        outputs[0] = outputArray;
    }

    void c2r(ArgumentList outputs, ArgumentList inputs) {
        auto handle = reinterpret_cast<PlanHandle*>(static_cast<uint64_t>(inputs[1][0]));
        TypedArray<std::complex<double>> inputArray = inputs[2];
        auto outputArray = ArrayFactory().createArray<double>(handle->realMatrixDims);
        fftw_execute_dft_c2r(handle->planInverse, reinterpret_cast<fftw_complex*>(inputArray.begin().operator->()), outputArray.begin().operator->());
        outputs[0] = outputArray;
    }

public:
    void operator()(ArgumentList outputs, ArgumentList inputs) override {
        if (inputs[0].getType() == matlab::data::ArrayType::CHAR) {
            matlab::data::CharArray commandArray = inputs[0];
            auto command = commandArray.toAscii();
            if (command == "create") createPlan(outputs, inputs);
            else if (command == "free") freePlan(outputs, inputs);
            else if (command == "r2c") r2c(outputs, inputs);
            else if (command == "c2r") c2r(outputs, inputs);
            else matlabPtr->feval(u"error", 0, {ArrayFactory().createScalar("Unknown command.")});
        } else {
            matlabPtr->feval(u"error", 0, {ArrayFactory().createScalar("Invalid input.")});
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
};
