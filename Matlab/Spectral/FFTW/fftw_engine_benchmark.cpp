#include "mex.hpp"
#include "mexAdapter.hpp"
#include <fftw3.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstdint>
#include <dlfcn.h>
#include <limits>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

using matlab::data::Array;
using matlab::data::ArrayDimensions;
using matlab::data::ArrayFactory;
using matlab::data::CharArray;
using matlab::data::TypedArray;
using matlab::mex::ArgumentList;

namespace {

using Clock = std::chrono::steady_clock;

struct ScratchBuffer {
    double* base = nullptr;
    double* data = nullptr;
};

struct PlanHandle {
    fftw_plan forwardReal = nullptr;
    fftw_plan forwardComplex = nullptr;
    fftw_plan inverseComplex = nullptr;
    fftw_plan inverseReal = nullptr;
    ArrayDimensions realDimensions;
    ArrayDimensions complexDimensions;
    std::vector<size_t> transformDimensions;
    std::string strategy;
    std::string alignmentMode;
    int forwardInputAlignment = 0;
    int forwardOutputAlignment = 0;
    int inverseInputAlignment = 0;
    int inverseOutputAlignment = 0;
    double planningSeconds = 0;
    bool planningLimitReached = false;
};

size_t product(const ArrayDimensions& dimensions) {
    return std::accumulate(dimensions.begin(), dimensions.end(), size_t{1}, std::multiplies<size_t>());
}

std::vector<size_t> strides(const ArrayDimensions& dimensions) {
    std::vector<size_t> values(dimensions.size(),1);
    for (size_t i = 1; i < dimensions.size(); ++i) {
        values[i] = values[i-1]*dimensions[i-1];
    }
    return values;
}

int checkedInt(size_t value, const char* label) {
    if (value > static_cast<size_t>(std::numeric_limits<int>::max())) {
        throw std::invalid_argument(std::string(label) + " exceeds the FFTW guru integer range.");
    }
    return static_cast<int>(value);
}

template <typename T>
const T* readPointer(const Array& input) {
    const TypedArray<T> typed = input;
    matlab::data::TypedIterator<const T> iterator(typed.begin());
    return iterator.operator->();
}

std::vector<size_t> numericVector(const Array& input) {
    TypedArray<double> values = input;
    std::vector<size_t> output;
    output.reserve(values.getNumberOfElements());
    for (double value : values) {
        if (!std::isfinite(value) || value < 1 || value != std::floor(value)) {
            throw std::invalid_argument("Dimension values must be positive integers.");
        }
        output.push_back(static_cast<size_t>(value));
    }
    return output;
}

std::string textValue(const Array& input) {
    CharArray value = input;
    return value.toAscii();
}

ArrayDimensions outputDimensions(const ArrayDimensions& realDimensions, const std::vector<size_t>& transformDimensions) {
    ArrayDimensions output = realDimensions;
    const size_t compressedDimension = transformDimensions.back();
    output[compressedDimension] = realDimensions[compressedDimension]/2 + 1;
    return output;
}

std::vector<fftw_iodim> transformIODimensions(const ArrayDimensions& realDimensions, const ArrayDimensions& complexDimensions, const std::vector<size_t>& transformDimensions, bool forward) {
    const std::vector<size_t> realStrides = strides(realDimensions);
    const std::vector<size_t> complexStrides = strides(complexDimensions);
    std::vector<fftw_iodim> dimensions;
    dimensions.reserve(transformDimensions.size());
    for (size_t dimension : transformDimensions) {
        fftw_iodim ioDimension;
        ioDimension.n = checkedInt(realDimensions[dimension],"Transform length");
        ioDimension.is = checkedInt(forward ? realStrides[dimension] : complexStrides[dimension],"Input stride");
        ioDimension.os = checkedInt(forward ? complexStrides[dimension] : realStrides[dimension],"Output stride");
        dimensions.push_back(ioDimension);
    }
    return dimensions;
}

std::vector<fftw_iodim> batchIODimensions(const ArrayDimensions& realDimensions, const ArrayDimensions& complexDimensions, const std::vector<size_t>& transformDimensions, bool forward) {
    const std::vector<size_t> realStrides = strides(realDimensions);
    const std::vector<size_t> complexStrides = strides(complexDimensions);
    std::vector<fftw_iodim> dimensions;
    for (size_t dimension = 0; dimension < realDimensions.size(); ++dimension) {
        if (std::find(transformDimensions.begin(),transformDimensions.end(),dimension) != transformDimensions.end()) {
            continue;
        }
        fftw_iodim ioDimension;
        ioDimension.n = checkedInt(realDimensions[dimension],"Batch length");
        ioDimension.is = checkedInt(forward ? realStrides[dimension] : complexStrides[dimension],"Batch input stride");
        ioDimension.os = checkedInt(forward ? complexStrides[dimension] : realStrides[dimension],"Batch output stride");
        dimensions.push_back(ioDimension);
    }
    return dimensions;
}

std::vector<fftw_iodim> stagedBatchDimensions(const ArrayDimensions& dimensions, size_t transformDimension) {
    const std::vector<size_t> dimensionStrides = strides(dimensions);
    std::vector<fftw_iodim> batches;
    for (size_t dimension = 0; dimension < dimensions.size(); ++dimension) {
        if (dimension == transformDimension) {
            continue;
        }
        batches.push_back(fftw_iodim{checkedInt(dimensions[dimension],"Batch length"),checkedInt(dimensionStrides[dimension],"Batch stride"),checkedInt(dimensionStrides[dimension],"Batch stride")});
    }
    return batches;
}

ScratchBuffer alignedBuffer(size_t nDoubles, int targetAlignment) {
    constexpr size_t padding = 64;
    ScratchBuffer buffer;
    buffer.base = fftw_alloc_real(nDoubles + padding);
    if (!buffer.base) {
        throw std::bad_alloc();
    }
    for (size_t offset = 0; offset < padding; ++offset) {
        if (fftw_alignment_of(buffer.base + offset) == targetAlignment) {
            buffer.data = buffer.base + offset;
            return buffer;
        }
    }
    fftw_free(buffer.base);
    throw std::runtime_error("Unable to manufacture the requested FFTW alignment class.");
}

void destroyPlan(PlanHandle* handle) {
    if (!handle) {
        return;
    }
    if (handle->forwardReal) fftw_destroy_plan(handle->forwardReal);
    if (handle->forwardComplex) fftw_destroy_plan(handle->forwardComplex);
    if (handle->inverseComplex) fftw_destroy_plan(handle->inverseComplex);
    if (handle->inverseReal) fftw_destroy_plan(handle->inverseReal);
    delete handle;
}

double elapsedSeconds(const Clock::time_point& start, const Clock::time_point& end) {
    return std::chrono::duration<double>(end-start).count();
}

bool alignmentClassesMatch(int expectedInput, int expectedOutput, int actualInput, int actualOutput) {
    return expectedInput == actualInput && expectedOutput == actualOutput;
}

} // namespace

class MexFunction : public matlab::mex::Function {
    ArrayFactory factory;
    std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();

    [[noreturn]] void fail(const std::string& identifier, const std::string& message) {
        matlabPtr->feval(u"error",0,{factory.createScalar(identifier),factory.createScalar(message)});
        throw std::runtime_error(message);
    }

    void require(bool condition, const std::string& identifier, const std::string& message) {
        if (!condition) {
            fail(identifier,message);
        }
    }

    PlanHandle* handleFrom(const Array& input) {
        const uint64_t value = static_cast<uint64_t>(input[0]);
        require(value != 0,"FFTWEngineBenchmark:InvalidPlan","The FFTW benchmark plan handle is invalid.");
        return reinterpret_cast<PlanHandle*>(value);
    }

    void validateDimensions(const ArrayDimensions& actual, const ArrayDimensions& expected, const std::string& label) {
        require(actual == expected,"FFTWEngineBenchmark:DimensionMismatch",label + " dimensions do not match the plan.");
    }

    void validateAlignment(const PlanHandle* handle, int inputAlignment, int outputAlignment, bool forward) {
        if (handle->alignmentMode == "unaligned") {
            return;
        }
        const int expectedInput = forward ? handle->forwardInputAlignment : handle->inverseInputAlignment;
        const int expectedOutput = forward ? handle->forwardOutputAlignment : handle->inverseOutputAlignment;
        require(alignmentClassesMatch(expectedInput,expectedOutput,inputAlignment,outputAlignment),"FFTWEngineBenchmark:AlignmentMismatch","The new-array alignment classes do not match the classes used to create this plan.");
    }

    fftw_plan createRankRealPlan(const PlanHandle& handle, double* realData, fftw_complex* complexData, unsigned flags, bool forward) {
        const auto transformDimensions = transformIODimensions(handle.realDimensions,handle.complexDimensions,handle.transformDimensions,forward);
        const auto batchDimensions = batchIODimensions(handle.realDimensions,handle.complexDimensions,handle.transformDimensions,forward);
        if (forward) {
            return fftw_plan_guru_dft_r2c(static_cast<int>(transformDimensions.size()),transformDimensions.data(),static_cast<int>(batchDimensions.size()),batchDimensions.data(),realData,complexData,flags);
        }
        return fftw_plan_guru_dft_c2r(static_cast<int>(transformDimensions.size()),transformDimensions.data(),static_cast<int>(batchDimensions.size()),batchDimensions.data(),complexData,realData,flags);
    }

    fftw_plan createStagedRealPlan(const PlanHandle& handle, double* realData, fftw_complex* complexData, unsigned flags, bool forward) {
        const size_t compressedDimension = handle.transformDimensions.back();
        const std::vector<size_t> realStrides = strides(handle.realDimensions);
        const std::vector<size_t> complexStrides = strides(handle.complexDimensions);
        fftw_iodim transformDimension{checkedInt(handle.realDimensions[compressedDimension],"Transform length"),checkedInt(forward ? realStrides[compressedDimension] : complexStrides[compressedDimension],"Input stride"),checkedInt(forward ? complexStrides[compressedDimension] : realStrides[compressedDimension],"Output stride")};
        const auto batches = batchIODimensions(handle.realDimensions,handle.complexDimensions,{compressedDimension},forward);
        if (forward) {
            return fftw_plan_guru_dft_r2c(1,&transformDimension,static_cast<int>(batches.size()),batches.data(),realData,complexData,flags);
        }
        return fftw_plan_guru_dft_c2r(1,&transformDimension,static_cast<int>(batches.size()),batches.data(),complexData,realData,flags);
    }

    fftw_plan createStagedComplexPlan(const PlanHandle& handle, fftw_complex* complexData, unsigned flags, bool forward) {
        const size_t complexDimension = handle.transformDimensions.front();
        const std::vector<size_t> complexStrides = strides(handle.complexDimensions);
        fftw_iodim transformDimension{checkedInt(handle.complexDimensions[complexDimension],"Transform length"),checkedInt(complexStrides[complexDimension],"Input stride"),checkedInt(complexStrides[complexDimension],"Output stride")};
        const auto batches = stagedBatchDimensions(handle.complexDimensions,complexDimension);
        return fftw_plan_guru_dft(1,&transformDimension,static_cast<int>(batches.size()),batches.data(),complexData,complexData,forward ? FFTW_FORWARD : FFTW_BACKWARD,flags);
    }

    void create(ArgumentList outputs, ArgumentList inputs) {
        require(inputs.size() == 10,"FFTWEngineBenchmark:InvalidInputCount","create expects dimensions, ordered transform dimensions, strategy, thread count, planner flags, alignment mode, time limit, real template, and spectrum template.");
        require(outputs.size() >= 1 && outputs.size() <= 8,"FFTWEngineBenchmark:InvalidOutputCount","create returns up to eight outputs.");

        const std::vector<size_t> dimensionValues = numericVector(inputs[1]);
        const std::vector<size_t> orderedTransformValues = numericVector(inputs[2]);
        require(orderedTransformValues.size() == 2,"FFTWEngineBenchmark:InvalidTransformDimensions","Exactly two ordered transform dimensions are required.");
        require(orderedTransformValues[0] != orderedTransformValues[1],"FFTWEngineBenchmark:RepeatedTransformDimension","Transform dimensions must be distinct.");

        auto handle = std::make_unique<PlanHandle>();
        handle->realDimensions.assign(dimensionValues.begin(),dimensionValues.end());
        handle->transformDimensions = orderedTransformValues;
        for (size_t& dimension : handle->transformDimensions) {
            require(dimension <= handle->realDimensions.size(),"FFTWEngineBenchmark:InvalidTransformDimension","A transform dimension exceeds the input rank.");
            --dimension;
        }
        handle->complexDimensions = outputDimensions(handle->realDimensions,handle->transformDimensions);
        handle->strategy = textValue(inputs[3]);
        handle->alignmentMode = textValue(inputs[6]);
        require(handle->strategy == "guru-rank2" || handle->strategy == "staged-r2c-c2c","FFTWEngineBenchmark:UnknownStrategy","Strategy must be guru-rank2 or staged-r2c-c2c.");
        require(handle->alignmentMode == "matched" || handle->alignmentMode == "unaligned","FFTWEngineBenchmark:UnknownAlignmentMode","Alignment mode must be matched or unaligned.");

        const int nThreads = static_cast<int>(static_cast<double>(inputs[4][0]));
        unsigned flags = static_cast<unsigned>(static_cast<double>(inputs[5][0]));
        const double timeLimit = static_cast<double>(inputs[7][0]);
        require(nThreads > 0,"FFTWEngineBenchmark:InvalidThreadCount","Thread count must be positive.");
        require(timeLimit > 0 && std::isfinite(timeLimit),"FFTWEngineBenchmark:InvalidTimeLimit","Planning time limit must be finite and positive.");
        if (handle->alignmentMode == "unaligned") {
            flags |= FFTW_UNALIGNED;
        }

        require(inputs[8].getType() == matlab::data::ArrayType::DOUBLE,"FFTWEngineBenchmark:InvalidRealTemplate","The real template must contain real doubles.");
        require(inputs[9].getType() == matlab::data::ArrayType::COMPLEX_DOUBLE,"FFTWEngineBenchmark:InvalidSpectrumTemplate","The spectrum template must contain complex doubles.");
        validateDimensions(inputs[8].getDimensions(),handle->realDimensions,"Real template");
        validateDimensions(inputs[9].getDimensions(),handle->complexDimensions,"Spectrum template");

        const double* realTemplate = readPointer<double>(inputs[8]);
        const std::complex<double>* spectrumTemplate = readPointer<std::complex<double>>(inputs[9]);
        handle->forwardInputAlignment = fftw_alignment_of(const_cast<double*>(realTemplate));
        handle->forwardOutputAlignment = fftw_alignment_of(reinterpret_cast<double*>(const_cast<std::complex<double>*>(spectrumTemplate)));
        handle->inverseInputAlignment = handle->forwardOutputAlignment;
        handle->inverseOutputAlignment = handle->forwardInputAlignment;

        ScratchBuffer realScratch;
        ScratchBuffer complexScratch;
        try {
            realScratch = alignedBuffer(product(handle->realDimensions),handle->forwardInputAlignment);
            complexScratch = alignedBuffer(2*product(handle->complexDimensions),handle->forwardOutputAlignment);
            fftw_complex* complexData = reinterpret_cast<fftw_complex*>(complexScratch.data);
            require(fftw_init_threads() != 0,"FFTWEngineBenchmark:ThreadInitializationFailed","FFTW thread initialization failed.");
            fftw_plan_with_nthreads(nThreads);
            fftw_set_timelimit(timeLimit);
            auto createTimedPlan = [&](auto&& createFunction) {
                const auto planStart = Clock::now();
                fftw_plan plan = createFunction();
                const double planSeconds = elapsedSeconds(planStart,Clock::now());
                handle->planningSeconds += planSeconds;
                handle->planningLimitReached = handle->planningLimitReached || planSeconds >= 0.95*timeLimit;
                return plan;
            };
            if (handle->strategy == "guru-rank2") {
                handle->forwardReal = createTimedPlan([&]() { return createRankRealPlan(*handle,realScratch.data,complexData,flags,true); });
                handle->inverseReal = createTimedPlan([&]() { return createRankRealPlan(*handle,realScratch.data,complexData,flags,false); });
            } else {
                handle->forwardReal = createTimedPlan([&]() { return createStagedRealPlan(*handle,realScratch.data,complexData,flags,true); });
                handle->forwardComplex = createTimedPlan([&]() { return createStagedComplexPlan(*handle,complexData,flags,true); });
                handle->inverseComplex = createTimedPlan([&]() { return createStagedComplexPlan(*handle,complexData,flags,false); });
                handle->inverseReal = createTimedPlan([&]() { return createStagedRealPlan(*handle,realScratch.data,complexData,flags,false); });
            }
            fftw_set_timelimit(FFTW_NO_TIMELIMIT);
            require(handle->forwardReal != nullptr && handle->inverseReal != nullptr,"FFTWEngineBenchmark:PlanCreationFailed","FFTW failed to create a real-transform plan.");
            if (handle->strategy == "staged-r2c-c2c") {
                require(handle->forwardComplex != nullptr && handle->inverseComplex != nullptr,"FFTWEngineBenchmark:PlanCreationFailed","FFTW failed to create a staged complex-transform plan.");
            }
        } catch (...) {
            fftw_set_timelimit(FFTW_NO_TIMELIMIT);
            if (realScratch.base) fftw_free(realScratch.base);
            if (complexScratch.base) fftw_free(complexScratch.base);
            destroyPlan(handle.release());
            throw;
        }
        fftw_free(realScratch.base);
        fftw_free(complexScratch.base);

        if (outputs.size() > 0) outputs[0] = factory.createScalar(reinterpret_cast<uint64_t>(handle.release()));
        PlanHandle* returnedHandle = reinterpret_cast<PlanHandle*>(static_cast<uint64_t>(outputs[0][0]));
        if (outputs.size() > 1) {
            auto dimensions = factory.createArray<double>({1,returnedHandle->complexDimensions.size()});
            std::copy(returnedHandle->complexDimensions.begin(),returnedHandle->complexDimensions.end(),dimensions.begin());
            outputs[1] = std::move(dimensions);
        }
        if (outputs.size() > 2) outputs[2] = factory.createScalar(1.0/static_cast<double>(returnedHandle->realDimensions[returnedHandle->transformDimensions[0]]*returnedHandle->realDimensions[returnedHandle->transformDimensions[1]]));
        if (outputs.size() > 3) outputs[3] = factory.createScalar(returnedHandle->planningSeconds);
        if (outputs.size() > 4) outputs[4] = factory.createScalar(returnedHandle->planningLimitReached);
        if (outputs.size() > 5) outputs[5] = factory.createScalar(static_cast<double>(returnedHandle->forwardInputAlignment));
        if (outputs.size() > 6) outputs[6] = factory.createScalar(static_cast<double>(returnedHandle->forwardOutputAlignment));
        if (outputs.size() > 7) outputs[7] = factory.createScalar(static_cast<double>(nThreads));
    }

    void free(ArgumentList inputs) {
        require(inputs.size() == 2,"FFTWEngineBenchmark:InvalidInputCount","free expects one plan handle.");
        destroyPlan(handleFrom(inputs[1]));
    }

    void r2c(ArgumentList outputs, ArgumentList inputs) {
        require(inputs.size() == 3 || inputs.size() == 4,"FFTWEngineBenchmark:InvalidInputCount","r2c expects a plan, real input, and optional preallocated spectrum.");
        require(outputs.size() >= 1 && outputs.size() <= 6,"FFTWEngineBenchmark:InvalidOutputCount","r2c returns the spectrum and timing/alignment metadata.");
        PlanHandle* handle = handleFrom(inputs[1]);
        require(inputs[2].getType() == matlab::data::ArrayType::DOUBLE,"FFTWEngineBenchmark:InvalidRealInput","r2c input must contain real doubles.");
        validateDimensions(inputs[2].getDimensions(),handle->realDimensions,"Real input");
        const double* inputPointer = readPointer<double>(inputs[2]);

        TypedArray<std::complex<double>> output = inputs.size() == 4 ? TypedArray<std::complex<double>>(std::move(inputs[3])) : factory.createArray<std::complex<double>>(handle->complexDimensions);
        validateDimensions(output.getDimensions(),handle->complexDimensions,"Spectrum output");
        std::complex<double>* outputPointer = &(*output.begin());
        const int inputAlignment = fftw_alignment_of(const_cast<double*>(inputPointer));
        const int outputAlignment = fftw_alignment_of(reinterpret_cast<double*>(outputPointer));
        validateAlignment(handle,inputAlignment,outputAlignment,true);

        const auto pipelineStart = Clock::now();
        const auto kernelStart = Clock::now();
        fftw_execute_dft_r2c(handle->forwardReal,const_cast<double*>(inputPointer),reinterpret_cast<fftw_complex*>(outputPointer));
        if (handle->forwardComplex) {
            fftw_execute_dft(handle->forwardComplex,reinterpret_cast<fftw_complex*>(outputPointer),reinterpret_cast<fftw_complex*>(outputPointer));
        }
        const auto kernelEnd = Clock::now();
        const auto pipelineEnd = Clock::now();

        outputs[0] = std::move(output);
        if (outputs.size() > 1) outputs[1] = factory.createScalar(elapsedSeconds(kernelStart,kernelEnd));
        if (outputs.size() > 2) outputs[2] = factory.createScalar(elapsedSeconds(pipelineStart,pipelineEnd));
        if (outputs.size() > 3) outputs[3] = factory.createScalar(static_cast<double>(inputAlignment));
        if (outputs.size() > 4) outputs[4] = factory.createScalar(static_cast<double>(outputAlignment));
        if (outputs.size() > 5) outputs[5] = factory.createScalar(handle->alignmentMode == "matched");
    }

    void c2r(ArgumentList outputs, ArgumentList inputs) {
        require(inputs.size() == 4,"FFTWEngineBenchmark:InvalidInputCount","c2r expects a plan, destructive spectrum input, and preallocated real output.");
        require(outputs.size() >= 2 && outputs.size() <= 7,"FFTWEngineBenchmark:InvalidOutputCount","c2r returns the destroyed spectrum, real output, and timing/alignment metadata.");
        PlanHandle* handle = handleFrom(inputs[1]);
        require(inputs[2].getType() == matlab::data::ArrayType::COMPLEX_DOUBLE,"FFTWEngineBenchmark:InvalidSpectrumInput","c2r input must contain complex doubles.");
        require(inputs[3].getType() == matlab::data::ArrayType::DOUBLE,"FFTWEngineBenchmark:InvalidRealOutput","c2r output must contain real doubles.");
        TypedArray<std::complex<double>> spectrum = std::move(inputs[2]);
        TypedArray<double> output = std::move(inputs[3]);
        validateDimensions(spectrum.getDimensions(),handle->complexDimensions,"Spectrum input");
        validateDimensions(output.getDimensions(),handle->realDimensions,"Real output");
        std::complex<double>* spectrumPointer = &(*spectrum.begin());
        double* outputPointer = &(*output.begin());
        const int inputAlignment = fftw_alignment_of(reinterpret_cast<double*>(spectrumPointer));
        const int outputAlignment = fftw_alignment_of(outputPointer);
        validateAlignment(handle,inputAlignment,outputAlignment,false);

        const auto pipelineStart = Clock::now();
        const auto kernelStart = Clock::now();
        if (handle->inverseComplex) {
            fftw_execute_dft(handle->inverseComplex,reinterpret_cast<fftw_complex*>(spectrumPointer),reinterpret_cast<fftw_complex*>(spectrumPointer));
        }
        fftw_execute_dft_c2r(handle->inverseReal,reinterpret_cast<fftw_complex*>(spectrumPointer),outputPointer);
        const auto kernelEnd = Clock::now();
        const auto pipelineEnd = Clock::now();

        outputs[0] = std::move(spectrum);
        outputs[1] = std::move(output);
        if (outputs.size() > 2) outputs[2] = factory.createScalar(elapsedSeconds(kernelStart,kernelEnd));
        if (outputs.size() > 3) outputs[3] = factory.createScalar(elapsedSeconds(pipelineStart,pipelineEnd));
        if (outputs.size() > 4) outputs[4] = factory.createScalar(static_cast<double>(inputAlignment));
        if (outputs.size() > 5) outputs[5] = factory.createScalar(static_cast<double>(outputAlignment));
        if (outputs.size() > 6) outputs[6] = factory.createScalar(handle->alignmentMode == "matched");
    }

    void forgetWisdom() {
        fftw_forget_wisdom();
    }

    void alignedCeiling(ArgumentList outputs, ArgumentList inputs) {
        require(inputs.size() == 4,"FFTWEngineBenchmark:InvalidInputCount","alignedCeiling expects a plan, real input, and sample count.");
        require(outputs.size() == 1,"FFTWEngineBenchmark:InvalidOutputCount","alignedCeiling returns raw FFTW-owned execution samples.");
        PlanHandle* handle = handleFrom(inputs[1]);
        require(inputs[2].getType() == matlab::data::ArrayType::DOUBLE,"FFTWEngineBenchmark:InvalidRealInput","alignedCeiling input must contain real doubles.");
        validateDimensions(inputs[2].getDimensions(),handle->realDimensions,"Real input");
        const int nSamples = static_cast<int>(static_cast<double>(inputs[3][0]));
        require(nSamples > 0,"FFTWEngineBenchmark:InvalidSampleCount","alignedCeiling sample count must be positive.");
        const double* inputPointer = readPointer<double>(inputs[2]);
        ScratchBuffer realScratch = alignedBuffer(product(handle->realDimensions),handle->forwardInputAlignment);
        ScratchBuffer complexScratch = alignedBuffer(2*product(handle->complexDimensions),handle->forwardOutputAlignment);
        fftw_complex* complexData = reinterpret_cast<fftw_complex*>(complexScratch.data);
        auto samples = factory.createArray<double>({1,static_cast<size_t>(nSamples)});
        auto sampleIterator = samples.begin();
        for (int iSample = 0; iSample < nSamples; ++iSample) {
            std::copy(inputPointer,inputPointer+product(handle->realDimensions),realScratch.data);
            const auto start = Clock::now();
            fftw_execute_dft_r2c(handle->forwardReal,realScratch.data,complexData);
            if (handle->forwardComplex) fftw_execute_dft(handle->forwardComplex,complexData,complexData);
            *sampleIterator++ = elapsedSeconds(start,Clock::now());
        }
        fftw_free(realScratch.base);
        fftw_free(complexScratch.base);
        outputs[0] = std::move(samples);
    }

    void alignmentSelfTest(ArgumentList outputs) {
        require(outputs.size() >= 2 && outputs.size() <= 3,"FFTWEngineBenchmark:InvalidOutputCount","alignmentSelfTest returns matchedAccepted, mismatchRejected, and optionally distinctPointerClassesObserved.");
        ScratchBuffer base;
        base.base = fftw_alloc_real(128);
        require(base.base != nullptr,"FFTWEngineBenchmark:AllocationFailed","Unable to allocate an FFTW alignment-test buffer.");
        base.data = base.base;
        const int matched = fftw_alignment_of(base.data);
        int mismatch = matched;
        for (size_t offset = 1; offset < 32 && mismatch == matched; ++offset) {
            mismatch = fftw_alignment_of(base.data + offset);
        }
        outputs[0] = factory.createScalar(alignmentClassesMatch(matched,matched,matched,matched));
        outputs[1] = factory.createScalar(!alignmentClassesMatch(matched,matched,matched+1,matched));
        if (outputs.size() > 2) outputs[2] = factory.createScalar(mismatch != matched);
        fftw_free(base.base);
    }

    void info(ArgumentList outputs) {
        require(outputs.size() == 2,"FFTWEngineBenchmark:InvalidOutputCount","info returns FFTW version and loaded library path.");
        Dl_info information{};
        std::string libraryPath = "unknown";
        if (dladdr(reinterpret_cast<void*>(fftw_execute),&information) != 0 && information.dli_fname) {
            libraryPath = information.dli_fname;
        }
        outputs[0] = factory.createScalar(std::string(fftw_version));
        outputs[1] = factory.createScalar(libraryPath);
    }

public:
    void operator()(ArgumentList outputs, ArgumentList inputs) override {
        require(!inputs.empty() && inputs[0].getType() == matlab::data::ArrayType::CHAR,"FFTWEngineBenchmark:InvalidCommand","The first input must be a command character vector.");
        const std::string command = textValue(inputs[0]);
        try {
            if (command == "create") create(outputs,inputs);
            else if (command == "free") free(inputs);
            else if (command == "r2c") r2c(outputs,inputs);
            else if (command == "c2r") c2r(outputs,inputs);
            else if (command == "forgetWisdom") forgetWisdom();
            else if (command == "alignedCeiling") alignedCeiling(outputs,inputs);
            else if (command == "alignmentSelfTest") alignmentSelfTest(outputs);
            else if (command == "info") info(outputs);
            else fail("FFTWEngineBenchmark:UnknownCommand","Unknown FFTW benchmark command: " + command);
        } catch (const matlab::engine::MATLABException&) {
            throw;
        } catch (const std::exception& exception) {
            fail("FFTWEngineBenchmark:NativeFailure",exception.what());
        }
    }
};
