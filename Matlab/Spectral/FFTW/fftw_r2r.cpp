#include "mex.hpp"
#include "mexAdapter.hpp"
#include "fftw_backend_support.hpp"
#include <fftw3.h>

#include <algorithm>
#include <complex>
#include <cstdint>
#include <dlfcn.h>
#include <memory>
#include <string>
#include <unordered_set>
#include <vector>

using matlab::data::Array;
using matlab::data::ArrayDimensions;
using matlab::data::ArrayFactory;
using matlab::data::TypedArray;
using matlab::mex::ArgumentList;

namespace {

using namespace fftw_backend;

struct Metrics {
    double allocationSeconds = 0;
    double wrapSeconds = 0;
    double kernelSeconds = 0;
    double normalizationSeconds = 0;
    double detachSeconds = 0;
    double internalSeconds = 0;
    double allocationCount = 0;
    double allocatedBytes = 0;
    double detectedCopyCount = 0;
    double detectedCopiedBytes = 0;
    double inputAlignment = 0;
    double outputAlignment = 0;
    uintptr_t inputPointer = 0;
    uintptr_t outputBefore = 0;
    uintptr_t outputMutable = 0;
    uintptr_t wrappedPointer = 0;
};

struct PlanHandle {
    fftw_plan forward = nullptr;
    fftw_plan inverse = nullptr;
    ArrayDimensions realDimensions;
    ArrayDimensions spectralDimensions;
    size_t transformDimension = 0;
    std::string transformType;
    std::string dataType;
    std::string alignmentMode;
    int forwardInputAlignment = 0;
    int forwardOutputAlignment = 0;
    int inverseInputAlignment = 0;
    int inverseOutputAlignment = 0;
    size_t realElements = 0;
    size_t spectralElements = 0;
    size_t scalarFactor = 1;
    double planningSeconds = 0;
    bool planningLimitReached = false;
    Metrics lastMetrics;
};

void destroyPlan(PlanHandle* handle) {
    if (!handle) return;
    if (handle->forward) fftw_destroy_plan(handle->forward);
    if (handle->inverse) fftw_destroy_plan(handle->inverse);
    delete handle;
}

} // namespace

class MexFunction : public matlab::mex::Function {
    ArrayFactory factory;
    std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();
    std::unordered_set<uint64_t> livePlans;
    size_t plansCreated = 0;
    size_t plansFreed = 0;
    size_t matlabBuffersCreated = 0;
    size_t matlabBuffersWrapped = 0;

    [[noreturn]] void fail(const std::string& identifier, const std::string& message) {
        matlabPtr->feval(u"error",0,{factory.createScalar(identifier),factory.createScalar(message)});
        throw std::runtime_error(message);
    }

    void require(bool condition, const std::string& identifier, const std::string& message) {
        if (!condition) fail(identifier,message);
    }

    PlanHandle* handleFrom(const Array& input) {
        const uint64_t token = static_cast<uint64_t>(input[0]);
        require(token != 0 && livePlans.count(token) == 1,"RealToRealTransform:InvalidPlan","The FFTW plan handle is invalid or has already been freed.");
        return reinterpret_cast<PlanHandle*>(token);
    }

    void validateDimensions(const ArrayDimensions& actual, const ArrayDimensions& expected, const std::string& label) {
        require(dimensionsEqualIgnoringTrailingSingletons(actual,expected),"RealToRealTransform:DimensionMismatch",label + " dimensions do not match the transform plan.");
    }

    void validateType(const Array& array, const PlanHandle* handle, const std::string& label) {
        const auto expected = handle->dataType == "real" ? matlab::data::ArrayType::DOUBLE : matlab::data::ArrayType::COMPLEX_DOUBLE;
        require(array.getType() == expected,"RealToRealTransform:DataTypeMismatch",label + " must match the configured dataType.");
    }

    const double* readDataPointer(const Array& array, const PlanHandle* handle) {
        if (handle->dataType == "real") return readPointer<double>(array);
        return reinterpret_cast<const double*>(readPointer<std::complex<double>>(array));
    }

    void validateAlignment(const PlanHandle* handle, int inputAlignment, int outputAlignment, bool forward) {
        if (handle->alignmentMode == "unaligned") return;
        const int expectedInput = forward ? handle->forwardInputAlignment : handle->inverseInputAlignment;
        const int expectedOutput = forward ? handle->forwardOutputAlignment : handle->inverseOutputAlignment;
        require(inputAlignment == expectedInput && outputAlignment == expectedOutput,"RealToRealTransform:AlignmentMismatch","The input or output alignment class differs from the matched FFTW plan; use alignmentMode=\"unaligned\" for arbitrary arrays.");
    }

    fftw_plan createPlan(const PlanHandle& handle, double* realData, double* spectralData, unsigned flags, bool forward) {
        const auto realStrides = strides(handle.realDimensions);
        const auto spectralStrides = strides(handle.spectralDimensions);
        const bool sine = handle.transformType == "sine";
        const size_t transformLength = sine ? handle.realDimensions[handle.transformDimension]-2 : handle.realDimensions[handle.transformDimension];
        const size_t inputStride = forward ? realStrides[handle.transformDimension] : spectralStrides[handle.transformDimension];
        const size_t outputStride = forward ? spectralStrides[handle.transformDimension] : realStrides[handle.transformDimension];
        fftw_iodim transform{checkedInt(transformLength,"Transform length"),checkedInt(handle.scalarFactor*inputStride,"Input stride"),checkedInt(handle.scalarFactor*outputStride,"Output stride")};

        std::vector<fftw_iodim> batches;
        if (handle.scalarFactor == 2) batches.push_back(fftw_iodim{2,1,1});
        for (size_t dimension = 0; dimension < handle.realDimensions.size(); ++dimension) {
            if (dimension == handle.transformDimension) continue;
            batches.push_back(fftw_iodim{
                checkedInt(handle.realDimensions[dimension],"Batch length"),
                checkedInt(handle.scalarFactor*(forward ? realStrides[dimension] : spectralStrides[dimension]),"Batch input stride"),
                checkedInt(handle.scalarFactor*(forward ? spectralStrides[dimension] : realStrides[dimension]),"Batch output stride")});
        }

        fftw_r2r_kind kind = sine ? FFTW_RODFT00 : FFTW_REDFT00;
        const size_t realOffset = sine ? handle.scalarFactor*realStrides[handle.transformDimension] : 0;
        double* input = forward ? realData+realOffset : spectralData;
        double* output = forward ? spectralData : realData+realOffset;
        return fftw_plan_guru_r2r(1,&transform,static_cast<int>(batches.size()),batches.data(),input,output,&kind,flags);
    }

    void create(ArgumentList outputs, ArgumentList inputs) {
        require(inputs.size() == 9,"RealToRealTransform:InvalidInputCount","create expects dimensions, transform dimension, transform type, data type, threads, planner flags, alignment mode, and planning time limit.");
        require(outputs.size() >= 1 && outputs.size() <= 7,"RealToRealTransform:InvalidOutputCount","create returns up to seven outputs.");
        auto handle = std::make_unique<PlanHandle>();
        const auto realDimensions = numericVector(inputs[1]);
        handle->realDimensions.assign(realDimensions.begin(),realDimensions.end());
        const auto transformDimensions = numericVector(inputs[2]);
        require(transformDimensions.size() == 1 && transformDimensions[0] <= handle->realDimensions.size(),"RealToRealTransform:InvalidTransformDimensions","Exactly one transform dimension within the input rank is required.");
        handle->transformDimension = transformDimensions[0]-1;
        require(handle->realDimensions[handle->transformDimension] > 2,"RealToRealTransform:TransformTooShort","DCT-I/DST-I transforms require at least three physical points.");
        handle->transformType = textValue(inputs[3]);
        handle->dataType = textValue(inputs[4]);
        require(handle->transformType == "cosine" || handle->transformType == "sine","RealToRealTransform:InvalidTransformType","Transform type must be cosine or sine.");
        require(handle->dataType == "real" || handle->dataType == "complex","RealToRealTransform:InvalidDataType","Data type must be real or complex.");
        handle->spectralDimensions = handle->realDimensions;
        if (handle->transformType == "sine") handle->spectralDimensions[handle->transformDimension] -= 2;
        handle->realElements = product(handle->realDimensions);
        handle->spectralElements = product(handle->spectralDimensions);
        handle->scalarFactor = handle->dataType == "complex" ? 2 : 1;

        const int nThreads = static_cast<int>(static_cast<double>(inputs[5][0]));
        unsigned flags = static_cast<unsigned>(static_cast<double>(inputs[6][0]));
        handle->alignmentMode = textValue(inputs[7]);
        const double timeLimit = static_cast<double>(inputs[8][0]);
        require(nThreads > 0,"RealToRealTransform:InvalidThreadCount","Thread count must be positive.");
        require(timeLimit > 0 && std::isfinite(timeLimit),"RealToRealTransform:InvalidTimeLimit","Planning time limit must be finite and positive.");
        require(handle->alignmentMode == "matched" || handle->alignmentMode == "unaligned","RealToRealTransform:InvalidAlignmentMode","Alignment mode must be matched or unaligned.");
        if (handle->alignmentMode == "unaligned") flags |= FFTW_UNALIGNED;

        const auto realStrides = strides(handle->realDimensions);
        const size_t realOffsetElements = handle->transformType == "sine" ? realStrides[handle->transformDimension] : 0;
        if (handle->dataType == "real") {
            auto realProbe = factory.createBuffer<double>(realOffsetElements+1);
            auto spectralProbe = factory.createBuffer<double>(1);
            handle->forwardInputAlignment = fftw_alignment_of(realProbe.get()+realOffsetElements);
            handle->forwardOutputAlignment = fftw_alignment_of(spectralProbe.get());
            handle->inverseInputAlignment = handle->forwardOutputAlignment;
            handle->inverseOutputAlignment = handle->forwardInputAlignment;
        } else {
            auto realProbe = factory.createBuffer<std::complex<double>>(realOffsetElements+1);
            auto spectralProbe = factory.createBuffer<std::complex<double>>(1);
            handle->forwardInputAlignment = fftw_alignment_of(reinterpret_cast<double*>(realProbe.get()+realOffsetElements));
            handle->forwardOutputAlignment = fftw_alignment_of(reinterpret_cast<double*>(spectralProbe.get()));
            handle->inverseInputAlignment = handle->forwardOutputAlignment;
            handle->inverseOutputAlignment = handle->forwardInputAlignment;
        }

        ScratchBuffer realScratch;
        ScratchBuffer spectralScratch;
        const size_t realOffsetScalars = handle->scalarFactor*realOffsetElements;
        try {
            realScratch = alignedBufferWithOffset(handle->scalarFactor*handle->realElements,realOffsetScalars,handle->forwardInputAlignment);
            spectralScratch = alignedBuffer(handle->scalarFactor*handle->spectralElements,handle->forwardOutputAlignment);
            require(fftw_init_threads() != 0,"RealToRealTransform:ThreadInitializationFailed","FFTW thread initialization failed.");
            fftw_plan_with_nthreads(nThreads);
            fftw_set_timelimit(timeLimit);
            auto createTimedPlan = [&](bool forward) {
                const auto start = Clock::now();
                fftw_plan plan = createPlan(*handle,realScratch.data,spectralScratch.data,flags,forward);
                const double elapsed = elapsedSeconds(start,Clock::now());
                handle->planningSeconds += elapsed;
                handle->planningLimitReached = handle->planningLimitReached || elapsed >= 0.95*timeLimit;
                return plan;
            };
            handle->forward = createTimedPlan(true);
            handle->inverse = createTimedPlan(false);
            fftw_set_timelimit(FFTW_NO_TIMELIMIT);
            require(handle->forward && handle->inverse,"RealToRealTransform:PlanCreationFailed","FFTW failed to create the forward or inverse plan.");
        } catch (...) {
            fftw_set_timelimit(FFTW_NO_TIMELIMIT);
            if (realScratch.base) fftw_free(realScratch.base);
            if (spectralScratch.base) fftw_free(spectralScratch.base);
            destroyPlan(handle.release());
            throw;
        }
        fftw_free(realScratch.base);
        fftw_free(spectralScratch.base);

        PlanHandle* raw = handle.release();
        const uint64_t token = reinterpret_cast<uint64_t>(raw);
        livePlans.insert(token);
        ++plansCreated;
        mexLock();
        try {
            outputs[0] = factory.createScalar(token);
            if (outputs.size() > 1) {
                auto dimensions = factory.createArray<double>({1,raw->spectralDimensions.size()});
                std::transform(raw->spectralDimensions.begin(),raw->spectralDimensions.end(),dimensions.begin(),[](size_t value) { return static_cast<double>(value); });
                outputs[1] = dimensions;
            }
            if (outputs.size() > 2) outputs[2] = factory.createScalar(1.0);
            if (outputs.size() > 3) outputs[3] = factory.createScalar(raw->planningSeconds);
            if (outputs.size() > 4) outputs[4] = factory.createScalar(raw->planningLimitReached);
            if (outputs.size() > 5) outputs[5] = factory.createScalar(static_cast<double>(raw->forwardInputAlignment));
            if (outputs.size() > 6) outputs[6] = factory.createScalar(static_cast<double>(raw->forwardOutputAlignment));
        } catch (...) {
            livePlans.erase(token);
            ++plansFreed;
            destroyPlan(raw);
            mexUnlock();
            throw;
        }
    }

    void freePlan(ArgumentList inputs) {
        require(inputs.size() == 2,"RealToRealTransform:InvalidInputCount","free expects a plan handle.");
        const uint64_t token = static_cast<uint64_t>(inputs[1][0]);
        if (token == 0 || livePlans.erase(token) == 0) return;
        destroyPlan(reinterpret_cast<PlanHandle*>(token));
        ++plansFreed;
        mexUnlock();
    }

    void normalize(PlanHandle* handle, const double* input, double* output, bool forward) {
        const auto realStrides = strides(handle->realDimensions);
        const auto spectralStrides = strides(handle->spectralDimensions);
        const size_t nPhysical = handle->realDimensions[handle->transformDimension];
        const size_t outputElements = forward ? handle->spectralElements : handle->realElements;
        const size_t outputStride = forward ? spectralStrides[handle->transformDimension] : realStrides[handle->transformDimension];
        const size_t components = handle->scalarFactor;
        if (handle->transformType == "cosine" && forward) {
            const double scale = 1.0/static_cast<double>(nPhysical-1);
            for (size_t element = 0; element < outputElements; ++element) {
                const size_t coordinate = (element/outputStride)%nPhysical;
                const double factor = coordinate == nPhysical-1 ? 0.5*scale : scale;
                for (size_t component = 0; component < components; ++component) output[components*element+component] *= factor;
            }
        } else if (handle->transformType == "cosine") {
            for (size_t element = 0; element < outputElements; ++element) {
                const size_t coordinate = (element/outputStride)%nPhysical;
                const size_t inputElement = element + (nPhysical-1-coordinate)*outputStride;
                const double sign = coordinate%2 == 0 ? 1.0 : -1.0;
                for (size_t component = 0; component < components; ++component) output[components*element+component] = 0.5*output[components*element+component] + 0.5*sign*input[components*inputElement+component];
            }
        } else if (forward) {
            const double scale = 1.0/static_cast<double>(nPhysical-1);
            for (size_t scalar = 0; scalar < components*outputElements; ++scalar) output[scalar] *= scale;
        } else {
            for (size_t element = 0; element < outputElements; ++element) {
                const size_t coordinate = (element/outputStride)%nPhysical;
                for (size_t component = 0; component < components; ++component) {
                    double& value = output[components*element+component];
                    value = coordinate == 0 || coordinate == nPhysical-1 ? 0.0 : 0.5*value;
                }
            }
        }
    }

    template <typename T>
    void executeTyped(ArgumentList outputs, ArgumentList inputs, PlanHandle* handle, bool forward, const ArrayDimensions& outputDimensions, size_t outputElements, const double* input, double* effectiveInput, const std::string& mode, Metrics& metrics) {
        const size_t byteCount = outputElements*sizeof(T);
        if (mode == "allocating") {
            const auto allocationStart = Clock::now();
            auto buffer = factory.createBuffer<T>(outputElements);
            ++matlabBuffersCreated;
            T* typedOutput = buffer.get();
            double* output = reinterpret_cast<double*>(typedOutput);
            metrics.allocationSeconds = elapsedSeconds(allocationStart,Clock::now());
            metrics.allocationCount = 1;
            metrics.allocatedBytes = static_cast<double>(byteCount);
            metrics.outputMutable = reinterpret_cast<uintptr_t>(typedOutput);
            const auto realStrides = strides(handle->realDimensions);
            const size_t outputOffset = !forward && handle->transformType == "sine" ? handle->scalarFactor*realStrides[handle->transformDimension] : 0;
            metrics.outputAlignment = fftw_alignment_of(output+outputOffset);
            validateAlignment(handle,static_cast<int>(metrics.inputAlignment),static_cast<int>(metrics.outputAlignment),forward);
            const auto kernelStart = Clock::now();
            fftw_execute_r2r(forward ? handle->forward : handle->inverse,effectiveInput,output+outputOffset);
            metrics.kernelSeconds = elapsedSeconds(kernelStart,Clock::now());
            const auto normalizationStart = Clock::now();
            normalize(handle,input,output,forward);
            metrics.normalizationSeconds = elapsedSeconds(normalizationStart,Clock::now());
            const auto wrapStart = Clock::now();
            auto result = factory.createArrayFromBuffer(outputDimensions,std::move(buffer));
            ++matlabBuffersWrapped;
            metrics.wrapSeconds = elapsedSeconds(wrapStart,Clock::now());
            metrics.wrappedPointer = metrics.outputMutable;
            outputs[0] = result;
        } else if (mode == "preallocated") {
            require(inputs.size() == 5,"RealToRealTransform:MissingOutput","Preallocated execution requires an output array.");
            validateType(inputs[4],handle,"Output");
            validateDimensions(inputs[4].getDimensions(),outputDimensions,"Output");
            metrics.outputBefore = reinterpret_cast<uintptr_t>(readPointer<T>(inputs[4]));
            const auto detachStart = Clock::now();
            TypedArray<T> result = std::move(inputs[4]);
            T* typedOutput = &(*result.begin());
            double* output = reinterpret_cast<double*>(typedOutput);
            metrics.detachSeconds = elapsedSeconds(detachStart,Clock::now());
            metrics.outputMutable = reinterpret_cast<uintptr_t>(typedOutput);
            metrics.wrappedPointer = metrics.outputMutable;
            if (metrics.outputBefore != metrics.outputMutable) {
                metrics.detectedCopyCount = 1;
                metrics.detectedCopiedBytes = static_cast<double>(byteCount);
                metrics.allocationCount = 1;
                metrics.allocatedBytes = static_cast<double>(byteCount);
            }
            const auto realStrides = strides(handle->realDimensions);
            const size_t outputOffset = !forward && handle->transformType == "sine" ? handle->scalarFactor*realStrides[handle->transformDimension] : 0;
            metrics.outputAlignment = fftw_alignment_of(output+outputOffset);
            validateAlignment(handle,static_cast<int>(metrics.inputAlignment),static_cast<int>(metrics.outputAlignment),forward);
            const auto kernelStart = Clock::now();
            fftw_execute_r2r(forward ? handle->forward : handle->inverse,effectiveInput,output+outputOffset);
            metrics.kernelSeconds = elapsedSeconds(kernelStart,Clock::now());
            const auto normalizationStart = Clock::now();
            normalize(handle,input,output,forward);
            metrics.normalizationSeconds = elapsedSeconds(normalizationStart,Clock::now());
            outputs[0] = result;
        } else {
            fail("RealToRealTransform:UnknownExecutionMode","Execution mode must be allocating or preallocated.");
        }
    }

    void execute(ArgumentList outputs, ArgumentList inputs, bool forward) {
        require(outputs.size() == 1,"RealToRealTransform:InvalidOutputCount","Transform execution returns one array.");
        require(inputs.size() == 4 || inputs.size() == 5,"RealToRealTransform:InvalidInputCount","Transform execution expects a plan, input, mode, and optional preallocated output.");
        PlanHandle* handle = handleFrom(inputs[1]);
        validateType(inputs[2],handle,"Input");
        const ArrayDimensions& inputDimensions = forward ? handle->realDimensions : handle->spectralDimensions;
        const ArrayDimensions& outputDimensions = forward ? handle->spectralDimensions : handle->realDimensions;
        validateDimensions(inputs[2].getDimensions(),inputDimensions,"Input");
        const double* input = readDataPointer(inputs[2],handle);
        const auto realStrides = strides(handle->realDimensions);
        const size_t inputOffset = forward && handle->transformType == "sine" ? handle->scalarFactor*realStrides[handle->transformDimension] : 0;
        double* effectiveInput = const_cast<double*>(input)+inputOffset;
        const std::string mode = textValue(inputs[3]);
        Metrics metrics;
        metrics.inputPointer = reinterpret_cast<uintptr_t>(input);
        metrics.inputAlignment = fftw_alignment_of(effectiveInput);
        const auto internalStart = Clock::now();
        const size_t outputElements = forward ? handle->spectralElements : handle->realElements;
        if (handle->dataType == "real") executeTyped<double>(outputs,inputs,handle,forward,outputDimensions,outputElements,input,effectiveInput,mode,metrics);
        else executeTyped<std::complex<double>>(outputs,inputs,handle,forward,outputDimensions,outputElements,input,effectiveInput,mode,metrics);
        metrics.internalSeconds = elapsedSeconds(internalStart,Clock::now());
        handle->lastMetrics = metrics;
    }

    void metrics(ArgumentList outputs, ArgumentList inputs) {
        require(outputs.size() >= 1 && outputs.size() <= 2,"RealToRealTransform:InvalidOutputCount","metrics returns values and optional pointer tokens.");
        PlanHandle* handle = handleFrom(inputs[1]);
        const Metrics& m = handle->lastMetrics;
        const std::vector<double> values = {m.allocationSeconds,m.wrapSeconds,m.kernelSeconds,m.normalizationSeconds,m.detachSeconds,m.internalSeconds,m.allocationCount,m.allocatedBytes,m.detectedCopyCount,m.detectedCopiedBytes,m.inputAlignment,m.outputAlignment,static_cast<double>(livePlans.size())};
        outputs[0] = factory.createArray<double>({1,values.size()},values.data(),values.data()+values.size());
        if (outputs.size() > 1) {
            const std::vector<uint64_t> pointers = {static_cast<uint64_t>(m.inputPointer),static_cast<uint64_t>(m.outputBefore),static_cast<uint64_t>(m.outputMutable),static_cast<uint64_t>(m.wrappedPointer)};
            outputs[1] = factory.createArray<uint64_t>({1,pointers.size()},pointers.data(),pointers.data()+pointers.size());
        }
    }

    void pointerToken(ArgumentList outputs, ArgumentList inputs) {
        require(outputs.size() == 1 && inputs.size() == 2,"RealToRealTransform:InvalidPointerCall","pointer expects one real or complex array.");
        uintptr_t pointer = 0;
        if (inputs[1].getType() == matlab::data::ArrayType::DOUBLE) pointer = reinterpret_cast<uintptr_t>(readPointer<double>(inputs[1]));
        else if (inputs[1].getType() == matlab::data::ArrayType::COMPLEX_DOUBLE) pointer = reinterpret_cast<uintptr_t>(readPointer<std::complex<double>>(inputs[1]));
        else fail("RealToRealTransform:InvalidPointerType","pointer accepts real or complex double arrays.");
        outputs[0] = factory.createScalar(static_cast<uint64_t>(pointer));
    }

    void lifetime(ArgumentList outputs) {
        require(outputs.size() == 1,"RealToRealTransform:InvalidOutputCount","lifetime returns one counter vector.");
        outputs[0] = factory.createArray<double>({1,6},{static_cast<double>(plansCreated),static_cast<double>(plansFreed),static_cast<double>(livePlans.size()),static_cast<double>(matlabBuffersCreated),static_cast<double>(matlabBuffersWrapped),static_cast<double>(livePlans.size())});
    }

    void resetLifetime(ArgumentList inputs) {
        require(inputs.size() == 1,"RealToRealTransform:InvalidInputCount","resetLifetime accepts no additional inputs.");
        require(livePlans.empty(),"RealToRealTransform:OutstandingPlans","Cannot reset lifetime counters while plans remain live.");
        plansCreated = 0;
        plansFreed = 0;
        matlabBuffersCreated = 0;
        matlabBuffersWrapped = 0;
    }

    void planInfo(ArgumentList outputs, ArgumentList inputs) {
        require(outputs.size() == 1,"RealToRealTransform:InvalidOutputCount","planInfo returns one numeric record.");
        PlanHandle* handle = handleFrom(inputs[1]);
        outputs[0] = factory.createArray<double>({1,9},{static_cast<double>(handle->realElements),static_cast<double>(handle->spectralElements),0.0,handle->planningSeconds,static_cast<double>(handle->planningLimitReached),static_cast<double>(handle->forwardInputAlignment),static_cast<double>(handle->forwardOutputAlignment),static_cast<double>(handle->inverseInputAlignment),static_cast<double>(handle->inverseOutputAlignment)});
    }

    void alignmentSelfTest(ArgumentList outputs) {
        require(outputs.size() >= 2 && outputs.size() <= 3,"RealToRealTransform:InvalidOutputCount","alignmentSelfTest returns matched acceptance, mismatch rejection, and optional unaligned acceptance.");
        double* base = fftw_alloc_real(128);
        require(base != nullptr,"RealToRealTransform:AllocationFailed","Unable to allocate alignment self-test storage.");
        const int matched = fftw_alignment_of(base);
        int mismatch = matched;
        for (size_t offset = 1; offset < 64 && mismatch == matched; ++offset) mismatch = fftw_alignment_of(base+offset);
        fftw_free(base);
        const int incompatible = mismatch == matched ? matched+1 : mismatch;
        outputs[0] = factory.createScalar(true);
        outputs[1] = factory.createScalar(incompatible != matched);
        if (outputs.size() > 2) outputs[2] = factory.createScalar(true);
    }

    void info(ArgumentList outputs) {
        Dl_info details{};
        const int found = dladdr(reinterpret_cast<void*>(reinterpret_cast<uintptr_t>(&fftw_execute)),&details);
        outputs[0] = factory.createScalar(fftw_version);
        if (outputs.size() > 1) outputs[1] = factory.createScalar(found && details.dli_fname ? details.dli_fname : "unknown");
    }

public:
    ~MexFunction() override {
        for (uint64_t token : livePlans) destroyPlan(reinterpret_cast<PlanHandle*>(token));
        livePlans.clear();
    }

    void operator()(ArgumentList outputs, ArgumentList inputs) override {
        require(!inputs.empty() && inputs[0].getType() == matlab::data::ArrayType::CHAR,"RealToRealTransform:InvalidCommand","The first input must be a command character vector.");
        const std::string command = textValue(inputs[0]);
        if (command == "create") create(outputs,inputs);
        else if (command == "free") freePlan(inputs);
        else if (command == "forward") execute(outputs,inputs,true);
        else if (command == "back") execute(outputs,inputs,false);
        else if (command == "metrics") metrics(outputs,inputs);
        else if (command == "pointer") pointerToken(outputs,inputs);
        else if (command == "lifetime") lifetime(outputs);
        else if (command == "resetLifetime") resetLifetime(inputs);
        else if (command == "planInfo") planInfo(outputs,inputs);
        else if (command == "alignmentSelfTest") alignmentSelfTest(outputs);
        else if (command == "forgetWisdom") fftw_forget_wisdom();
        else if (command == "info") info(outputs);
        else if (command == "noop") return;
        else fail("RealToRealTransform:UnknownCommand","Unknown FFTW r2r backend command.");
    }
};
