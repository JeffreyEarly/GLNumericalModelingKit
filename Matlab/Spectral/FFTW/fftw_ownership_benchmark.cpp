#include "mex.hpp"
#include "mexAdapter.hpp"
#include <fftw3.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstdint>
#include <cstring>
#include <dlfcn.h>
#include <limits>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

#ifdef __APPLE__
#include <mach/mach.h>
#include <malloc/malloc.h>
#endif

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

struct Metrics {
    double allocationSeconds = 0;
    double wrapSeconds = 0;
    double memcpySeconds = 0;
    double detachSeconds = 0;
    double kernelSeconds = 0;
    double internalSeconds = 0;
    double allocationCount = 0;
    double allocatedBytes = 0;
    double explicitCopyCount = 0;
    double explicitCopiedBytes = 0;
    double detectedCopyCount = 0;
    double detectedCopiedBytes = 0;
    double inputAlignment = 0;
    double outputAlignment = 0;
    double destroyedInput = 0;
    double persistentScratchBytes = 0;
    uintptr_t inputBefore = 0;
    uintptr_t inputMutable = 0;
    uintptr_t outputBefore = 0;
    uintptr_t outputMutable = 0;
    uintptr_t wrappedPointer = 0;
};

struct PlanHandle {
    fftw_plan forward = nullptr;
    fftw_plan inverse = nullptr;
    ArrayDimensions realDimensions;
    ArrayDimensions complexDimensions;
    std::vector<size_t> transformDimensions;
    std::string alignmentMode;
    int forwardInputAlignment = 0;
    int forwardOutputAlignment = 0;
    int inverseInputAlignment = 0;
    int inverseOutputAlignment = 0;
    double* complexScratchBase = nullptr;
    fftw_complex* complexScratch = nullptr;
    size_t realElements = 0;
    size_t complexElements = 0;
    double planningSeconds = 0;
    bool planningLimitReached = false;
    Metrics lastMetrics;
};

size_t product(const ArrayDimensions& dimensions) {
    return std::accumulate(dimensions.begin(),dimensions.end(),size_t{1},std::multiplies<size_t>());
}

std::vector<size_t> strides(const ArrayDimensions& dimensions) {
    std::vector<size_t> values(dimensions.size(),1);
    for (size_t i = 1; i < dimensions.size(); ++i) values[i] = values[i-1]*dimensions[i-1];
    return values;
}

int checkedInt(size_t value, const char* label) {
    if (value > static_cast<size_t>(std::numeric_limits<int>::max())) throw std::invalid_argument(std::string(label) + " exceeds the FFTW guru integer range.");
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
        if (!std::isfinite(value) || value < 1 || value != std::floor(value)) throw std::invalid_argument("Dimension values must be positive integers.");
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
    const auto realStrides = strides(realDimensions);
    const auto complexStrides = strides(complexDimensions);
    std::vector<fftw_iodim> dimensions;
    for (size_t dimension : transformDimensions) {
        dimensions.push_back(fftw_iodim{checkedInt(realDimensions[dimension],"Transform length"),checkedInt(forward ? realStrides[dimension] : complexStrides[dimension],"Input stride"),checkedInt(forward ? complexStrides[dimension] : realStrides[dimension],"Output stride")});
    }
    return dimensions;
}

std::vector<fftw_iodim> batchIODimensions(const ArrayDimensions& realDimensions, const ArrayDimensions& complexDimensions, const std::vector<size_t>& transformDimensions, bool forward) {
    const auto realStrides = strides(realDimensions);
    const auto complexStrides = strides(complexDimensions);
    std::vector<fftw_iodim> dimensions;
    for (size_t dimension = 0; dimension < realDimensions.size(); ++dimension) {
        if (std::find(transformDimensions.begin(),transformDimensions.end(),dimension) != transformDimensions.end()) continue;
        dimensions.push_back(fftw_iodim{checkedInt(realDimensions[dimension],"Batch length"),checkedInt(forward ? realStrides[dimension] : complexStrides[dimension],"Batch input stride"),checkedInt(forward ? complexStrides[dimension] : realStrides[dimension],"Batch output stride")});
    }
    return dimensions;
}

ScratchBuffer alignedBuffer(size_t nDoubles, int targetAlignment) {
    constexpr size_t padding = 64;
    ScratchBuffer buffer;
    buffer.base = fftw_alloc_real(nDoubles + padding);
    if (!buffer.base) throw std::bad_alloc();
    for (size_t offset = 0; offset < padding; ++offset) {
        if (fftw_alignment_of(buffer.base + offset) == targetAlignment) {
            buffer.data = buffer.base + offset;
            return buffer;
        }
    }
    fftw_free(buffer.base);
    throw std::runtime_error("Unable to manufacture the requested FFTW alignment class.");
}

double elapsedSeconds(const Clock::time_point& start, const Clock::time_point& end) {
    return std::chrono::duration<double>(end-start).count();
}

void destroyPlan(PlanHandle* handle) {
    if (!handle) return;
    if (handle->forward) fftw_destroy_plan(handle->forward);
    if (handle->inverse) fftw_destroy_plan(handle->inverse);
    if (handle->complexScratchBase) fftw_free(handle->complexScratchBase);
    delete handle;
}

} // namespace

class MexFunction : public matlab::mex::Function {
    ArrayFactory factory;
    std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();
    static MexFunction* instance;
    static size_t fftwBuffersCreated;
    static size_t fftwBuffersFreed;
    static size_t fftwBuffersOutstanding;

    static void deleteFftwBuffer(std::complex<double>* pointer) {
        fftw_free(pointer);
        ++fftwBuffersFreed;
        if (fftwBuffersOutstanding > 0) --fftwBuffersOutstanding;
        if (instance) instance->mexUnlock();
    }

    [[noreturn]] void fail(const std::string& identifier, const std::string& message) {
        matlabPtr->feval(u"error",0,{factory.createScalar(identifier),factory.createScalar(message)});
        throw std::runtime_error(message);
    }

    void require(bool condition, const std::string& identifier, const std::string& message) {
        if (!condition) fail(identifier,message);
    }

    PlanHandle* handleFrom(const Array& input) {
        const uint64_t value = static_cast<uint64_t>(input[0]);
        require(value != 0,"FFTWMexOwnership:InvalidPlan","The ownership benchmark plan handle is invalid.");
        return reinterpret_cast<PlanHandle*>(value);
    }

    void validateDimensions(const ArrayDimensions& actual, const ArrayDimensions& expected, const std::string& label) {
        require(actual == expected,"FFTWMexOwnership:DimensionMismatch",label + " dimensions do not match the plan.");
    }

    void validateAlignment(const PlanHandle* handle, int inputAlignment, int outputAlignment, bool forward) {
        if (handle->alignmentMode == "unaligned") return;
        const int expectedInput = forward ? handle->forwardInputAlignment : handle->inverseInputAlignment;
        const int expectedOutput = forward ? handle->forwardOutputAlignment : handle->inverseOutputAlignment;
        require(inputAlignment == expectedInput && outputAlignment == expectedOutput,"FFTWMexOwnership:AlignmentMismatch","The new-array alignment classes do not match the plan.");
    }

    fftw_plan createPlan(const PlanHandle& handle, double* realData, fftw_complex* complexData, unsigned flags, bool forward) {
        const auto transforms = transformIODimensions(handle.realDimensions,handle.complexDimensions,handle.transformDimensions,forward);
        const auto batches = batchIODimensions(handle.realDimensions,handle.complexDimensions,handle.transformDimensions,forward);
        if (forward) return fftw_plan_guru_dft_r2c(static_cast<int>(transforms.size()),transforms.data(),static_cast<int>(batches.size()),batches.data(),realData,complexData,flags);
        return fftw_plan_guru_dft_c2r(static_cast<int>(transforms.size()),transforms.data(),static_cast<int>(batches.size()),batches.data(),complexData,realData,flags);
    }

    void create(ArgumentList outputs, ArgumentList inputs) {
        require(inputs.size() == 8 || inputs.size() == 9,"FFTWMexOwnership:InvalidInputCount","create expects dimensions, transform order, threads, planner flags, alignment mode, real and spectrum templates, and an optional planning time limit.");
        auto handle = std::make_unique<PlanHandle>();
        const auto realDimensions = numericVector(inputs[1]);
        handle->realDimensions.assign(realDimensions.begin(),realDimensions.end());
        const auto ordered = numericVector(inputs[2]);
        require(ordered.size() == 2 && ordered[0] != ordered[1],"FFTWMexOwnership:InvalidTransformOrder","Exactly two distinct transform dimensions are required.");
        handle->transformDimensions = ordered;
        for (size_t& dimension : handle->transformDimensions) {
            require(dimension <= handle->realDimensions.size(),"FFTWMexOwnership:InvalidTransformDimension","A transform dimension exceeds the input rank.");
            --dimension;
        }
        handle->complexDimensions = outputDimensions(handle->realDimensions,handle->transformDimensions);
        handle->realElements = product(handle->realDimensions);
        handle->complexElements = product(handle->complexDimensions);
        const int nThreads = static_cast<int>(static_cast<double>(inputs[3][0]));
        unsigned flags = static_cast<unsigned>(static_cast<double>(inputs[4][0]));
        const double timeLimit = inputs.size() == 9 ? static_cast<double>(inputs[8][0]) : 10.0;
        require(nThreads > 0,"FFTWMexOwnership:InvalidThreadCount","Thread count must be positive.");
        require(timeLimit > 0 && std::isfinite(timeLimit),"FFTWMexOwnership:InvalidTimeLimit","Planning time limit must be finite and positive.");
        handle->alignmentMode = textValue(inputs[5]);
        require(handle->alignmentMode == "matched" || handle->alignmentMode == "unaligned","FFTWMexOwnership:InvalidAlignmentMode","Alignment mode must be matched or unaligned.");
        if (handle->alignmentMode == "unaligned") flags |= FFTW_UNALIGNED;
        validateDimensions(inputs[6].getDimensions(),handle->realDimensions,"Real template");
        validateDimensions(inputs[7].getDimensions(),handle->complexDimensions,"Spectrum template");
        const double* realTemplate = readPointer<double>(inputs[6]);
        const std::complex<double>* spectrumTemplate = readPointer<std::complex<double>>(inputs[7]);
        handle->forwardInputAlignment = fftw_alignment_of(const_cast<double*>(realTemplate));
        handle->forwardOutputAlignment = fftw_alignment_of(reinterpret_cast<double*>(const_cast<std::complex<double>*>(spectrumTemplate)));
        handle->inverseInputAlignment = handle->forwardOutputAlignment;
        handle->inverseOutputAlignment = handle->forwardInputAlignment;
        ScratchBuffer realScratch;
        ScratchBuffer complexScratch;
        try {
            realScratch = alignedBuffer(handle->realElements,handle->forwardInputAlignment);
            complexScratch = alignedBuffer(2*handle->complexElements,handle->forwardOutputAlignment);
            require(fftw_init_threads() != 0,"FFTWMexOwnership:ThreadInitializationFailed","FFTW thread initialization failed.");
            fftw_plan_with_nthreads(nThreads);
            fftw_set_timelimit(timeLimit);
            auto createTimedPlan = [&](bool forward) {
                const auto start = Clock::now();
                fftw_plan plan = createPlan(*handle,realScratch.data,reinterpret_cast<fftw_complex*>(complexScratch.data),flags,forward);
                const double elapsed = elapsedSeconds(start,Clock::now());
                handle->planningSeconds += elapsed;
                handle->planningLimitReached = handle->planningLimitReached || elapsed >= 0.95*timeLimit;
                return plan;
            };
            handle->forward = createTimedPlan(true);
            handle->inverse = createTimedPlan(false);
            fftw_set_timelimit(FFTW_NO_TIMELIMIT);
            require(handle->forward && handle->inverse,"FFTWMexOwnership:PlanCreationFailed","FFTW failed to create the ownership benchmark plans.");
            handle->complexScratchBase = complexScratch.base;
            handle->complexScratch = reinterpret_cast<fftw_complex*>(complexScratch.data);
            complexScratch.base = nullptr;
        } catch (...) {
            fftw_set_timelimit(FFTW_NO_TIMELIMIT);
            if (realScratch.base) fftw_free(realScratch.base);
            if (complexScratch.base) fftw_free(complexScratch.base);
            destroyPlan(handle.release());
            throw;
        }
        fftw_free(realScratch.base);
        PlanHandle* raw = handle.release();
        outputs[0] = factory.createScalar(reinterpret_cast<uint64_t>(raw));
        auto dimensions = factory.createArray<double>({1,raw->complexDimensions.size()});
        std::transform(raw->complexDimensions.begin(),raw->complexDimensions.end(),dimensions.begin(),[](size_t value) { return static_cast<double>(value); });
        if (outputs.size() > 1) outputs[1] = dimensions;
        double scale = 1;
        for (size_t dimension : raw->transformDimensions) scale /= static_cast<double>(raw->realDimensions[dimension]);
        if (outputs.size() > 2) outputs[2] = factory.createScalar(scale);
        if (outputs.size() > 3) outputs[3] = factory.createScalar(raw->planningSeconds);
        if (outputs.size() > 4) outputs[4] = factory.createScalar(static_cast<double>(raw->forwardInputAlignment));
        if (outputs.size() > 5) outputs[5] = factory.createScalar(static_cast<double>(raw->forwardOutputAlignment));
        if (outputs.size() > 6) outputs[6] = factory.createScalar(raw->planningLimitReached);
    }

    void freePlan(ArgumentList, ArgumentList inputs) {
        require(inputs.size() == 2,"FFTWMexOwnership:InvalidInputCount","free expects a plan handle.");
        destroyPlan(handleFrom(inputs[1]));
    }

    void forward(ArgumentList outputs, ArgumentList inputs) {
        require(inputs.size() == 4 || inputs.size() == 5,"FFTWMexOwnership:InvalidInputCount","forward expects a plan, real input, mode, and optional caller output.");
        PlanHandle* handle = handleFrom(inputs[1]);
        validateDimensions(inputs[2].getDimensions(),handle->realDimensions,"Forward input");
        const std::string mode = textValue(inputs[3]);
        const double* input = readPointer<double>(inputs[2]);
        Metrics metrics;
        metrics.inputBefore = reinterpret_cast<uintptr_t>(input);
        metrics.inputMutable = metrics.inputBefore;
        metrics.inputAlignment = fftw_alignment_of(const_cast<double*>(input));
        metrics.persistentScratchBytes = static_cast<double>(handle->complexElements*sizeof(fftw_complex));
        const auto internalStart = Clock::now();

        if (mode == "factory-array") {
            const auto start = Clock::now();
            auto output = factory.createArray<std::complex<double>>(handle->complexDimensions);
            std::complex<double>* pointer = &(*output.begin());
            metrics.allocationSeconds = elapsedSeconds(start,Clock::now());
            metrics.allocationCount = 1;
            metrics.allocatedBytes = handle->complexElements*sizeof(std::complex<double>);
            metrics.outputMutable = reinterpret_cast<uintptr_t>(pointer);
            metrics.wrappedPointer = metrics.outputMutable;
            metrics.outputAlignment = fftw_alignment_of(reinterpret_cast<double*>(pointer));
            validateAlignment(handle,static_cast<int>(metrics.inputAlignment),static_cast<int>(metrics.outputAlignment),true);
            const auto kernelStart = Clock::now();
            fftw_execute_dft_r2c(handle->forward,const_cast<double*>(input),reinterpret_cast<fftw_complex*>(pointer));
            metrics.kernelSeconds = elapsedSeconds(kernelStart,Clock::now());
            outputs[0] = output;
        } else if (mode == "caller-direct") {
            require(inputs.size() == 5,"FFTWMexOwnership:MissingCallerOutput","caller-direct requires a preallocated output.");
            validateDimensions(inputs[4].getDimensions(),handle->complexDimensions,"Caller output");
            metrics.outputBefore = reinterpret_cast<uintptr_t>(readPointer<std::complex<double>>(inputs[4]));
            const auto detachStart = Clock::now();
            TypedArray<std::complex<double>> output = std::move(inputs[4]);
            std::complex<double>* pointer = &(*output.begin());
            metrics.detachSeconds = elapsedSeconds(detachStart,Clock::now());
            metrics.outputMutable = reinterpret_cast<uintptr_t>(pointer);
            metrics.wrappedPointer = metrics.outputMutable;
            if (metrics.outputBefore != metrics.outputMutable) {
                metrics.detectedCopyCount = 1;
                metrics.detectedCopiedBytes = handle->complexElements*sizeof(std::complex<double>);
                metrics.allocationCount = 1;
                metrics.allocatedBytes = metrics.detectedCopiedBytes;
            }
            metrics.outputAlignment = fftw_alignment_of(reinterpret_cast<double*>(pointer));
            validateAlignment(handle,static_cast<int>(metrics.inputAlignment),static_cast<int>(metrics.outputAlignment),true);
            const auto kernelStart = Clock::now();
            fftw_execute_dft_r2c(handle->forward,const_cast<double*>(input),reinterpret_cast<fftw_complex*>(pointer));
            metrics.kernelSeconds = elapsedSeconds(kernelStart,Clock::now());
            outputs[0] = output;
        } else if (mode == "matlab-buffer") {
            const auto allocationStart = Clock::now();
            auto buffer = factory.createBuffer<std::complex<double>>(handle->complexElements);
            std::complex<double>* pointer = buffer.get();
            metrics.allocationSeconds = elapsedSeconds(allocationStart,Clock::now());
            metrics.allocationCount = 1;
            metrics.allocatedBytes = handle->complexElements*sizeof(std::complex<double>);
            metrics.outputMutable = reinterpret_cast<uintptr_t>(pointer);
            metrics.outputAlignment = fftw_alignment_of(reinterpret_cast<double*>(pointer));
            validateAlignment(handle,static_cast<int>(metrics.inputAlignment),static_cast<int>(metrics.outputAlignment),true);
            const auto kernelStart = Clock::now();
            fftw_execute_dft_r2c(handle->forward,const_cast<double*>(input),reinterpret_cast<fftw_complex*>(pointer));
            metrics.kernelSeconds = elapsedSeconds(kernelStart,Clock::now());
            const auto wrapStart = Clock::now();
            auto output = factory.createArrayFromBuffer(handle->complexDimensions,std::move(buffer));
            metrics.wrapSeconds = elapsedSeconds(wrapStart,Clock::now());
            metrics.wrappedPointer = metrics.outputMutable;
            outputs[0] = output;
        } else if (mode == "fftw-buffer") {
            static_assert(sizeof(std::complex<double>) == sizeof(fftw_complex),"FFTW and std::complex<double> storage sizes differ.");
            const auto allocationStart = Clock::now();
            auto* pointer = reinterpret_cast<std::complex<double>*>(fftw_alloc_complex(handle->complexElements));
            metrics.allocationSeconds = elapsedSeconds(allocationStart,Clock::now());
            require(pointer != nullptr,"FFTWMexOwnership:AllocationFailed","fftw_alloc_complex failed.");
            metrics.allocationCount = 1;
            metrics.allocatedBytes = handle->complexElements*sizeof(std::complex<double>);
            metrics.outputMutable = reinterpret_cast<uintptr_t>(pointer);
            metrics.outputAlignment = fftw_alignment_of(reinterpret_cast<double*>(pointer));
            try {
                validateAlignment(handle,static_cast<int>(metrics.inputAlignment),static_cast<int>(metrics.outputAlignment),true);
                const auto kernelStart = Clock::now();
                fftw_execute_dft_r2c(handle->forward,const_cast<double*>(input),reinterpret_cast<fftw_complex*>(pointer));
                metrics.kernelSeconds = elapsedSeconds(kernelStart,Clock::now());
            } catch (...) {
                fftw_free(pointer);
                throw;
            }
            const auto wrapStart = Clock::now();
            matlab::data::buffer_ptr_t<std::complex<double>> buffer(pointer,&MexFunction::deleteFftwBuffer);
            mexLock();
            ++fftwBuffersCreated;
            ++fftwBuffersOutstanding;
            try {
                auto output = factory.createArrayFromBuffer(handle->complexDimensions,std::move(buffer));
                metrics.wrapSeconds = elapsedSeconds(wrapStart,Clock::now());
                metrics.wrappedPointer = metrics.outputMutable;
                outputs[0] = output;
            } catch (...) {
                if (buffer) {
                    buffer.reset();
                }
                throw;
            }
        } else {
            fail("FFTWMexOwnership:UnknownForwardMode","Unknown forward ownership mode.");
        }
        metrics.internalSeconds = elapsedSeconds(internalStart,Clock::now());
        handle->lastMetrics = metrics;
    }

    void inversePreserving(ArgumentList outputs, ArgumentList inputs) {
        require(inputs.size() == 4 || inputs.size() == 5,"FFTWMexOwnership:InvalidInputCount","inversePreserving expects a plan, spectrum, mode, and optional real output.");
        PlanHandle* handle = handleFrom(inputs[1]);
        validateDimensions(inputs[2].getDimensions(),handle->complexDimensions,"Preserving inverse input");
        const std::string mode = textValue(inputs[3]);
        const std::complex<double>* input = readPointer<std::complex<double>>(inputs[2]);
        Metrics metrics;
        metrics.inputBefore = reinterpret_cast<uintptr_t>(input);
        metrics.inputMutable = metrics.inputBefore;
        metrics.inputAlignment = fftw_alignment_of(reinterpret_cast<double*>(const_cast<std::complex<double>*>(input)));
        metrics.persistentScratchBytes = static_cast<double>(handle->complexElements*sizeof(fftw_complex));
        const auto internalStart = Clock::now();
        const auto copyStart = Clock::now();
        std::memcpy(handle->complexScratch,input,handle->complexElements*sizeof(fftw_complex));
        metrics.memcpySeconds = elapsedSeconds(copyStart,Clock::now());
        metrics.explicitCopyCount = 1;
        metrics.explicitCopiedBytes = handle->complexElements*sizeof(fftw_complex);
        double* outputPointer = nullptr;
        std::unique_ptr<TypedArray<double>> output;
        if (mode == "allocating") {
            const auto allocationStart = Clock::now();
            output = std::make_unique<TypedArray<double>>(factory.createArray<double>(handle->realDimensions));
            outputPointer = &(*output->begin());
            metrics.allocationSeconds = elapsedSeconds(allocationStart,Clock::now());
            metrics.allocationCount = 1;
            metrics.allocatedBytes = handle->realElements*sizeof(double);
        } else if (mode == "preallocated") {
            require(inputs.size() == 5,"FFTWMexOwnership:MissingCallerOutput","preallocated preserving inverse requires a real output.");
            validateDimensions(inputs[4].getDimensions(),handle->realDimensions,"Preserving inverse output");
            metrics.outputBefore = reinterpret_cast<uintptr_t>(readPointer<double>(inputs[4]));
            const auto detachStart = Clock::now();
            output = std::make_unique<TypedArray<double>>(std::move(inputs[4]));
            outputPointer = &(*output->begin());
            metrics.detachSeconds = elapsedSeconds(detachStart,Clock::now());
            metrics.outputMutable = reinterpret_cast<uintptr_t>(outputPointer);
            if (metrics.outputBefore != metrics.outputMutable) {
                metrics.detectedCopyCount = 1;
                metrics.detectedCopiedBytes = handle->realElements*sizeof(double);
                metrics.allocationCount = 1;
                metrics.allocatedBytes = metrics.detectedCopiedBytes;
            }
        } else {
            fail("FFTWMexOwnership:UnknownInverseMode","Unknown preserving inverse mode.");
        }
        metrics.outputMutable = reinterpret_cast<uintptr_t>(outputPointer);
        metrics.wrappedPointer = metrics.outputMutable;
        metrics.outputAlignment = fftw_alignment_of(outputPointer);
        validateAlignment(handle,fftw_alignment_of(reinterpret_cast<double*>(handle->complexScratch)),static_cast<int>(metrics.outputAlignment),false);
        const auto kernelStart = Clock::now();
        fftw_execute_dft_c2r(handle->inverse,handle->complexScratch,outputPointer);
        metrics.kernelSeconds = elapsedSeconds(kernelStart,Clock::now());
        outputs[0] = *output;
        metrics.internalSeconds = elapsedSeconds(internalStart,Clock::now());
        handle->lastMetrics = metrics;
    }

    void inverseDestructive(ArgumentList outputs, ArgumentList inputs) {
        require(inputs.size() == 4,"FFTWMexOwnership:InvalidInputCount","inverseDestructive expects a plan, spectrum, and real output.");
        PlanHandle* handle = handleFrom(inputs[1]);
        validateDimensions(inputs[2].getDimensions(),handle->complexDimensions,"Destructive inverse input");
        validateDimensions(inputs[3].getDimensions(),handle->realDimensions,"Destructive inverse output");
        Metrics metrics;
        metrics.persistentScratchBytes = static_cast<double>(handle->complexElements*sizeof(fftw_complex));
        metrics.inputBefore = reinterpret_cast<uintptr_t>(readPointer<std::complex<double>>(inputs[2]));
        metrics.outputBefore = reinterpret_cast<uintptr_t>(readPointer<double>(inputs[3]));
        const auto internalStart = Clock::now();
        const auto detachStart = Clock::now();
        TypedArray<std::complex<double>> spectrum = std::move(inputs[2]);
        TypedArray<double> output = std::move(inputs[3]);
        std::complex<double>* spectrumPointer = &(*spectrum.begin());
        double* outputPointer = &(*output.begin());
        metrics.detachSeconds = elapsedSeconds(detachStart,Clock::now());
        metrics.inputMutable = reinterpret_cast<uintptr_t>(spectrumPointer);
        metrics.outputMutable = reinterpret_cast<uintptr_t>(outputPointer);
        metrics.wrappedPointer = metrics.outputMutable;
        if (metrics.inputBefore != metrics.inputMutable) {
            metrics.detectedCopyCount += 1;
            metrics.detectedCopiedBytes += handle->complexElements*sizeof(std::complex<double>);
            metrics.allocationCount += 1;
            metrics.allocatedBytes += handle->complexElements*sizeof(std::complex<double>);
        }
        if (metrics.outputBefore != metrics.outputMutable) {
            metrics.detectedCopyCount += 1;
            metrics.detectedCopiedBytes += handle->realElements*sizeof(double);
            metrics.allocationCount += 1;
            metrics.allocatedBytes += handle->realElements*sizeof(double);
        }
        metrics.inputAlignment = fftw_alignment_of(reinterpret_cast<double*>(spectrumPointer));
        metrics.outputAlignment = fftw_alignment_of(outputPointer);
        validateAlignment(handle,static_cast<int>(metrics.inputAlignment),static_cast<int>(metrics.outputAlignment),false);
        const auto kernelStart = Clock::now();
        fftw_execute_dft_c2r(handle->inverse,reinterpret_cast<fftw_complex*>(spectrumPointer),outputPointer);
        metrics.kernelSeconds = elapsedSeconds(kernelStart,Clock::now());
        metrics.destroyedInput = 1;
        outputs[0] = spectrum;
        outputs[1] = output;
        metrics.internalSeconds = elapsedSeconds(internalStart,Clock::now());
        handle->lastMetrics = metrics;
    }

    void metrics(ArgumentList outputs, ArgumentList inputs) {
        PlanHandle* handle = handleFrom(inputs[1]);
        const Metrics& m = handle->lastMetrics;
        const std::vector<double> values = {m.allocationSeconds,m.wrapSeconds,m.memcpySeconds,m.detachSeconds,m.kernelSeconds,m.internalSeconds,m.allocationCount,m.allocatedBytes,m.explicitCopyCount,m.explicitCopiedBytes,m.detectedCopyCount,m.detectedCopiedBytes,m.inputAlignment,m.outputAlignment,m.destroyedInput,m.persistentScratchBytes,static_cast<double>(fftwBuffersOutstanding)};
        outputs[0] = factory.createArray<double>({1,values.size()},values.data(),values.data()+values.size());
        if (outputs.size() > 1) {
            const std::vector<uint64_t> pointers = {static_cast<uint64_t>(m.inputBefore),static_cast<uint64_t>(m.inputMutable),static_cast<uint64_t>(m.outputBefore),static_cast<uint64_t>(m.outputMutable),static_cast<uint64_t>(m.wrappedPointer)};
            outputs[1] = factory.createArray<uint64_t>({1,pointers.size()},pointers.data(),pointers.data()+pointers.size());
        }
    }

    void pointerToken(ArgumentList outputs, ArgumentList inputs) {
        require(inputs.size() == 2,"FFTWMexOwnership:InvalidInputCount","pointer expects one array.");
        uintptr_t pointer = 0;
        if (inputs[1].getType() == matlab::data::ArrayType::DOUBLE) pointer = reinterpret_cast<uintptr_t>(readPointer<double>(inputs[1]));
        else if (inputs[1].getType() == matlab::data::ArrayType::COMPLEX_DOUBLE) pointer = reinterpret_cast<uintptr_t>(readPointer<std::complex<double>>(inputs[1]));
        else fail("FFTWMexOwnership:InvalidPointerType","pointer accepts real or complex double arrays.");
        outputs[0] = factory.createScalar(static_cast<uint64_t>(pointer));
    }

    void lifetime(ArgumentList outputs) {
        outputs[0] = factory.createArray<double>({1,3},{static_cast<double>(fftwBuffersCreated),static_cast<double>(fftwBuffersFreed),static_cast<double>(fftwBuffersOutstanding)});
    }

    void resetLifetime(ArgumentList) {
        require(fftwBuffersOutstanding == 0,"FFTWMexOwnership:OutstandingBuffers","Cannot reset lifetime counters while FFTW-owned arrays remain alive.");
        fftwBuffersCreated = 0;
        fftwBuffersFreed = 0;
    }

    void info(ArgumentList outputs) {
        Dl_info details{};
        const int found = dladdr(reinterpret_cast<void*>(reinterpret_cast<uintptr_t>(&fftw_execute)),&details);
        outputs[0] = factory.createScalar(fftw_version);
        if (outputs.size() > 1) outputs[1] = factory.createScalar(found && details.dli_fname ? details.dli_fname : "unknown");
    }

    void processMemory(ArgumentList outputs) {
        double allocatorBytes = std::numeric_limits<double>::quiet_NaN();
        double residentBytes = std::numeric_limits<double>::quiet_NaN();
#ifdef __APPLE__
        malloc_statistics_t statistics{};
        malloc_zone_statistics(malloc_default_zone(),&statistics);
        allocatorBytes = static_cast<double>(statistics.size_in_use);
        mach_task_basic_info_data_t info{};
        mach_msg_type_number_t count = MACH_TASK_BASIC_INFO_COUNT;
        if (task_info(mach_task_self(),MACH_TASK_BASIC_INFO,reinterpret_cast<task_info_t>(&info),&count) == KERN_SUCCESS) residentBytes = static_cast<double>(info.resident_size);
#endif
        outputs[0] = factory.createArray<double>({1,2},{allocatorBytes,residentBytes});
    }

    void alignmentSelfTest(ArgumentList outputs) {
        double* base = fftw_alloc_real(128);
        require(base != nullptr,"FFTWMexOwnership:AllocationFailed","Unable to allocate alignment self-test storage.");
        const int matched = fftw_alignment_of(base);
        fftw_free(base);
        outputs[0] = factory.createScalar(true);
        if (outputs.size() > 1) outputs[1] = factory.createScalar(matched != matched+1);
        if (outputs.size() > 2) outputs[2] = factory.createScalar(true);
    }

    void forgetWisdom(ArgumentList inputs) {
        require(inputs.size() == 1,"FFTWMexOwnership:InvalidInputCount","forgetWisdom accepts no additional inputs.");
        fftw_forget_wisdom();
    }

public:
    MexFunction() { instance = this; }
    ~MexFunction() override { if (instance == this) instance = nullptr; }

    void operator()(ArgumentList outputs, ArgumentList inputs) override {
        require(!inputs.empty() && inputs[0].getType() == matlab::data::ArrayType::CHAR,"FFTWMexOwnership:InvalidCommand","The first input must be a command character vector.");
        const std::string command = textValue(inputs[0]);
        if (command == "create") create(outputs,inputs);
        else if (command == "free") freePlan(outputs,inputs);
        else if (command == "forward") forward(outputs,inputs);
        else if (command == "inversePreserving") inversePreserving(outputs,inputs);
        else if (command == "inverseDestructive") inverseDestructive(outputs,inputs);
        else if (command == "metrics") metrics(outputs,inputs);
        else if (command == "pointer") pointerToken(outputs,inputs);
        else if (command == "lifetime") lifetime(outputs);
        else if (command == "resetLifetime") resetLifetime(inputs);
        else if (command == "info") info(outputs);
        else if (command == "memory") processMemory(outputs);
        else if (command == "alignmentSelfTest") alignmentSelfTest(outputs);
        else if (command == "forgetWisdom") forgetWisdom(inputs);
        else if (command == "noop") return;
        else fail("FFTWMexOwnership:UnknownCommand","Unknown ownership benchmark command.");
    }
};

MexFunction* MexFunction::instance = nullptr;
size_t MexFunction::fftwBuffersCreated = 0;
size_t MexFunction::fftwBuffersFreed = 0;
size_t MexFunction::fftwBuffersOutstanding = 0;
