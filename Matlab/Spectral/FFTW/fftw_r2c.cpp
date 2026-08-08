#include "mex.hpp"
#include "mexAdapter.hpp"
#include "fftw_backend_support.hpp"
#include <fftw3.h>

#include <algorithm>
#include <complex>
#include <cstdint>
#include <cstring>
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
        const uint64_t value = static_cast<uint64_t>(input[0]);
        require(value != 0 && livePlans.count(value) == 1,"RealToComplexTransform:InvalidPlan","The FFTW plan handle is invalid or has already been freed.");
        return reinterpret_cast<PlanHandle*>(value);
    }

    void validateDimensions(const ArrayDimensions& actual, const ArrayDimensions& expected, const std::string& label) {
        require(dimensionsEqualIgnoringTrailingSingletons(actual,expected),"RealToComplexTransform:DimensionMismatch",label + " dimensions do not match the transform plan.");
    }

    void validateAlignment(const PlanHandle* handle, int inputAlignment, int outputAlignment, bool forward) {
        if (handle->alignmentMode == "unaligned") return;
        const int expectedInput = forward ? handle->forwardInputAlignment : handle->inverseInputAlignment;
        const int expectedOutput = forward ? handle->forwardOutputAlignment : handle->inverseOutputAlignment;
        require(inputAlignment == expectedInput && outputAlignment == expectedOutput,"RealToComplexTransform:AlignmentMismatch","The input or output alignment class differs from the matched FFTW plan; rebuild with alignmentMode=\"unaligned\" for arbitrary arrays.");
    }

    fftw_plan createPlan(const PlanHandle& handle, double* realData, fftw_complex* complexData, unsigned flags, bool forward) {
        const auto transforms = transformIODimensions(handle.realDimensions,handle.complexDimensions,handle.transformDimensions,forward);
        const auto batches = batchIODimensions(handle.realDimensions,handle.complexDimensions,handle.transformDimensions,forward);
        if (forward) return fftw_plan_guru_dft_r2c(static_cast<int>(transforms.size()),transforms.data(),static_cast<int>(batches.size()),batches.data(),realData,complexData,flags);
        return fftw_plan_guru_dft_c2r(static_cast<int>(transforms.size()),transforms.data(),static_cast<int>(batches.size()),batches.data(),complexData,realData,flags);
    }

    void create(ArgumentList outputs, ArgumentList inputs) {
        require(inputs.size() == 7,"RealToComplexTransform:InvalidInputCount","create expects dimensions, ordered transform dimensions, threads, planner flags, alignment mode, and planning time limit.");
        require(outputs.size() >= 1 && outputs.size() <= 7,"RealToComplexTransform:InvalidOutputCount","create returns up to seven outputs.");
        auto handle = std::make_unique<PlanHandle>();
        const auto realDimensions = numericVector(inputs[1]);
        handle->realDimensions.assign(realDimensions.begin(),realDimensions.end());
        auto ordered = numericVector(inputs[2]);
        require(!ordered.empty() && hasDistinctDimensions(ordered),"RealToComplexTransform:InvalidTransformDimensions","Transform dimensions must be a nonempty ordered list of distinct dimensions.");
        handle->transformDimensions = ordered;
        for (size_t& dimension : handle->transformDimensions) {
            require(dimension <= handle->realDimensions.size(),"RealToComplexTransform:InvalidTransformDimension","A transform dimension exceeds the input rank.");
            --dimension;
            require(handle->realDimensions[dimension] > 1,"RealToComplexTransform:SingletonTransformDimension","Transform dimensions must have length greater than one.");
        }
        handle->complexDimensions = outputDimensions(handle->realDimensions,handle->transformDimensions);
        handle->realElements = product(handle->realDimensions);
        handle->complexElements = product(handle->complexDimensions);
        const int nThreads = static_cast<int>(static_cast<double>(inputs[3][0]));
        unsigned flags = static_cast<unsigned>(static_cast<double>(inputs[4][0]));
        handle->alignmentMode = textValue(inputs[5]);
        const double timeLimit = static_cast<double>(inputs[6][0]);
        require(nThreads > 0,"RealToComplexTransform:InvalidThreadCount","Thread count must be positive.");
        require(timeLimit > 0 && std::isfinite(timeLimit),"RealToComplexTransform:InvalidTimeLimit","Planning time limit must be finite and positive.");
        require(handle->alignmentMode == "matched" || handle->alignmentMode == "unaligned","RealToComplexTransform:InvalidAlignmentMode","Alignment mode must be matched or unaligned.");
        if (handle->alignmentMode == "unaligned") flags |= FFTW_UNALIGNED;

        auto realProbe = factory.createArray<double>({1,1});
        auto complexProbe = factory.createBuffer<std::complex<double>>(1);
        handle->forwardInputAlignment = fftw_alignment_of(&(*realProbe.begin()));
        handle->forwardOutputAlignment = fftw_alignment_of(reinterpret_cast<double*>(complexProbe.get()));
        handle->inverseInputAlignment = handle->forwardOutputAlignment;
        handle->inverseOutputAlignment = handle->forwardInputAlignment;

        ScratchBuffer realScratch;
        ScratchBuffer complexScratch;
        try {
            realScratch = alignedBuffer(handle->realElements,handle->forwardInputAlignment);
            complexScratch = alignedBuffer(2*handle->complexElements,handle->forwardOutputAlignment);
            require(fftw_init_threads() != 0,"RealToComplexTransform:ThreadInitializationFailed","FFTW thread initialization failed.");
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
            require(handle->forward && handle->inverse,"RealToComplexTransform:PlanCreationFailed","FFTW failed to create the forward or inverse plan.");
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
        const uint64_t token = reinterpret_cast<uint64_t>(raw);
        livePlans.insert(token);
        ++plansCreated;
        mexLock();
        try {
            outputs[0] = factory.createScalar(token);
            if (outputs.size() > 1) {
                auto dimensions = factory.createArray<double>({1,raw->complexDimensions.size()});
                std::transform(raw->complexDimensions.begin(),raw->complexDimensions.end(),dimensions.begin(),[](size_t value) { return static_cast<double>(value); });
                outputs[1] = dimensions;
            }
            double scale = 1;
            for (size_t dimension : raw->transformDimensions) scale /= static_cast<double>(raw->realDimensions[dimension]);
            if (outputs.size() > 2) outputs[2] = factory.createScalar(scale);
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
        require(inputs.size() == 2,"RealToComplexTransform:InvalidInputCount","free expects a plan handle.");
        const uint64_t token = static_cast<uint64_t>(inputs[1][0]);
        if (token == 0 || livePlans.erase(token) == 0) return;
        destroyPlan(reinterpret_cast<PlanHandle*>(token));
        ++plansFreed;
        mexUnlock();
    }

    void forward(ArgumentList outputs, ArgumentList inputs) {
        require(outputs.size() == 1,"RealToComplexTransform:InvalidOutputCount","forward returns one half spectrum.");
        require(inputs.size() == 4 || inputs.size() == 5,"RealToComplexTransform:InvalidInputCount","forward expects a plan, real input, mode, and optional preallocated output.");
        PlanHandle* handle = handleFrom(inputs[1]);
        require(inputs[2].getType() == matlab::data::ArrayType::DOUBLE,"RealToComplexTransform:InvalidRealInput","Forward input must be a real double array.");
        validateDimensions(inputs[2].getDimensions(),handle->realDimensions,"Forward input");
        const double* input = readPointer<double>(inputs[2]);
        const std::string mode = textValue(inputs[3]);
        Metrics metrics;
        metrics.inputBefore = reinterpret_cast<uintptr_t>(input);
        metrics.inputMutable = metrics.inputBefore;
        metrics.inputAlignment = fftw_alignment_of(const_cast<double*>(input));
        metrics.persistentScratchBytes = static_cast<double>(handle->complexElements*sizeof(fftw_complex));
        const auto internalStart = Clock::now();

        if (mode == "allocating") {
            const auto allocationStart = Clock::now();
            auto buffer = factory.createBuffer<std::complex<double>>(handle->complexElements);
            ++matlabBuffersCreated;
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
            ++matlabBuffersWrapped;
            metrics.wrapSeconds = elapsedSeconds(wrapStart,Clock::now());
            metrics.wrappedPointer = metrics.outputMutable;
            outputs[0] = output;
        } else if (mode == "preallocated") {
            require(inputs.size() == 5,"RealToComplexTransform:MissingOutput","Preallocated forward execution requires an output array.");
            require(inputs[4].getType() == matlab::data::ArrayType::COMPLEX_DOUBLE,"RealToComplexTransform:InvalidComplexOutput","Forward output must be a complex double array.");
            validateDimensions(inputs[4].getDimensions(),handle->complexDimensions,"Forward output");
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
        } else {
            fail("RealToComplexTransform:UnknownForwardMode","Forward mode must be allocating or preallocated.");
        }
        metrics.internalSeconds = elapsedSeconds(internalStart,Clock::now());
        handle->lastMetrics = metrics;
    }

    void inversePreserving(ArgumentList outputs, ArgumentList inputs) {
        require(outputs.size() == 1,"RealToComplexTransform:InvalidOutputCount","A preserving inverse returns one real array.");
        require(inputs.size() == 4 || inputs.size() == 5,"RealToComplexTransform:InvalidInputCount","inversePreserving expects a plan, spectrum, mode, and optional preallocated output.");
        PlanHandle* handle = handleFrom(inputs[1]);
        require(inputs[2].getType() == matlab::data::ArrayType::COMPLEX_DOUBLE,"RealToComplexTransform:InvalidComplexInput","Inverse input must be a complex double array.");
        validateDimensions(inputs[2].getDimensions(),handle->complexDimensions,"Preserving inverse input");
        const std::complex<double>* input = readPointer<std::complex<double>>(inputs[2]);
        const std::string mode = textValue(inputs[3]);
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
            require(inputs.size() == 5,"RealToComplexTransform:MissingOutput","Preallocated inverse execution requires a real output array.");
            require(inputs[4].getType() == matlab::data::ArrayType::DOUBLE,"RealToComplexTransform:InvalidRealOutput","Inverse output must be a real double array.");
            validateDimensions(inputs[4].getDimensions(),handle->realDimensions,"Preserving inverse output");
            metrics.outputBefore = reinterpret_cast<uintptr_t>(readPointer<double>(inputs[4]));
            const auto detachStart = Clock::now();
            output = std::make_unique<TypedArray<double>>(std::move(inputs[4]));
            outputPointer = &(*output->begin());
            metrics.detachSeconds = elapsedSeconds(detachStart,Clock::now());
            if (metrics.outputBefore != reinterpret_cast<uintptr_t>(outputPointer)) {
                metrics.detectedCopyCount = 1;
                metrics.detectedCopiedBytes = handle->realElements*sizeof(double);
                metrics.allocationCount = 1;
                metrics.allocatedBytes = metrics.detectedCopiedBytes;
            }
        } else {
            fail("RealToComplexTransform:UnknownInverseMode","Inverse mode must be allocating or preallocated.");
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
        require(outputs.size() == 2,"RealToComplexTransform:InvalidOutputCount","A destructive inverse returns the destroyed spectrum and real output.");
        require(inputs.size() == 4,"RealToComplexTransform:InvalidInputCount","inverseDestructive expects a plan, spectrum, and preallocated real output.");
        PlanHandle* handle = handleFrom(inputs[1]);
        require(inputs[2].getType() == matlab::data::ArrayType::COMPLEX_DOUBLE,"RealToComplexTransform:InvalidComplexInput","Inverse input must be a complex double array.");
        require(inputs[3].getType() == matlab::data::ArrayType::DOUBLE,"RealToComplexTransform:InvalidRealOutput","Inverse output must be a real double array.");
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
        require(outputs.size() >= 1 && outputs.size() <= 2,"RealToComplexTransform:InvalidOutputCount","metrics returns values and optional pointer tokens.");
        PlanHandle* handle = handleFrom(inputs[1]);
        const Metrics& m = handle->lastMetrics;
        const std::vector<double> values = {m.allocationSeconds,m.wrapSeconds,m.memcpySeconds,m.detachSeconds,m.kernelSeconds,m.internalSeconds,m.allocationCount,m.allocatedBytes,m.explicitCopyCount,m.explicitCopiedBytes,m.detectedCopyCount,m.detectedCopiedBytes,m.inputAlignment,m.outputAlignment,m.destroyedInput,m.persistentScratchBytes,static_cast<double>(livePlans.size())};
        outputs[0] = factory.createArray<double>({1,values.size()},values.data(),values.data()+values.size());
        if (outputs.size() > 1) {
            const std::vector<uint64_t> pointers = {static_cast<uint64_t>(m.inputBefore),static_cast<uint64_t>(m.inputMutable),static_cast<uint64_t>(m.outputBefore),static_cast<uint64_t>(m.outputMutable),static_cast<uint64_t>(m.wrappedPointer)};
            outputs[1] = factory.createArray<uint64_t>({1,pointers.size()},pointers.data(),pointers.data()+pointers.size());
        }
    }

    void pointerToken(ArgumentList outputs, ArgumentList inputs) {
        require(outputs.size() == 1 && inputs.size() == 2,"RealToComplexTransform:InvalidPointerCall","pointer expects one array and returns one token.");
        uintptr_t pointer = 0;
        if (inputs[1].getType() == matlab::data::ArrayType::DOUBLE) pointer = reinterpret_cast<uintptr_t>(readPointer<double>(inputs[1]));
        else if (inputs[1].getType() == matlab::data::ArrayType::COMPLEX_DOUBLE) pointer = reinterpret_cast<uintptr_t>(readPointer<std::complex<double>>(inputs[1]));
        else fail("RealToComplexTransform:InvalidPointerType","pointer accepts real or complex double arrays.");
        outputs[0] = factory.createScalar(static_cast<uint64_t>(pointer));
    }

    void lifetime(ArgumentList outputs) {
        require(outputs.size() == 1,"RealToComplexTransform:InvalidOutputCount","lifetime returns one counter vector.");
        outputs[0] = factory.createArray<double>({1,6},{static_cast<double>(plansCreated),static_cast<double>(plansFreed),static_cast<double>(livePlans.size()),static_cast<double>(matlabBuffersCreated),static_cast<double>(matlabBuffersWrapped),static_cast<double>(livePlans.size())});
    }

    void resetLifetime(ArgumentList inputs) {
        require(inputs.size() == 1,"RealToComplexTransform:InvalidInputCount","resetLifetime accepts no additional inputs.");
        require(livePlans.empty(),"RealToComplexTransform:OutstandingPlans","Cannot reset lifetime counters while plans remain live.");
        plansCreated = 0;
        plansFreed = 0;
        matlabBuffersCreated = 0;
        matlabBuffersWrapped = 0;
    }

    void planInfo(ArgumentList outputs, ArgumentList inputs) {
        require(outputs.size() == 1,"RealToComplexTransform:InvalidOutputCount","planInfo returns one numeric record.");
        PlanHandle* handle = handleFrom(inputs[1]);
        outputs[0] = factory.createArray<double>({1,7},{static_cast<double>(handle->realElements),static_cast<double>(handle->complexElements),static_cast<double>(handle->complexElements*sizeof(fftw_complex)),handle->planningSeconds,static_cast<double>(handle->planningLimitReached),static_cast<double>(handle->forwardInputAlignment),static_cast<double>(handle->forwardOutputAlignment)});
    }

    void alignmentSelfTest(ArgumentList outputs) {
        require(outputs.size() >= 2 && outputs.size() <= 3,"RealToComplexTransform:InvalidOutputCount","alignmentSelfTest returns matched acceptance, mismatch rejection, and optional unaligned acceptance.");
        double* base = fftw_alloc_real(128);
        require(base != nullptr,"RealToComplexTransform:AllocationFailed","Unable to allocate alignment self-test storage.");
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
        require(!inputs.empty() && inputs[0].getType() == matlab::data::ArrayType::CHAR,"RealToComplexTransform:InvalidCommand","The first input must be a command character vector.");
        const std::string command = textValue(inputs[0]);
        if (command == "create") create(outputs,inputs);
        else if (command == "free") freePlan(inputs);
        else if (command == "forward") forward(outputs,inputs);
        else if (command == "inversePreserving") inversePreserving(outputs,inputs);
        else if (command == "inverseDestructive") inverseDestructive(outputs,inputs);
        else if (command == "metrics") metrics(outputs,inputs);
        else if (command == "pointer") pointerToken(outputs,inputs);
        else if (command == "lifetime") lifetime(outputs);
        else if (command == "resetLifetime") resetLifetime(inputs);
        else if (command == "planInfo") planInfo(outputs,inputs);
        else if (command == "alignmentSelfTest") alignmentSelfTest(outputs);
        else if (command == "forgetWisdom") fftw_forget_wisdom();
        else if (command == "info") info(outputs);
        else if (command == "noop") return;
        else fail("RealToComplexTransform:UnknownCommand","Unknown FFTW r2c backend command.");
    }
};
