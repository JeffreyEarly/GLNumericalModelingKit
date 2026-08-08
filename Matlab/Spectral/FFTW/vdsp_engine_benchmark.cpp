#include "mex.hpp"
#include "mexAdapter.hpp"

#include <Accelerate/Accelerate.h>
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

struct PlanHandle {
    ArrayDimensions realDimensions;
    ArrayDimensions complexDimensions;
    std::vector<size_t> transformDimensions;
    std::string strategy;
    vDSP_DFT_SetupD realForward = nullptr;
    vDSP_DFT_SetupD realInverse = nullptr;
    vDSP_DFT_SetupD complexForward = nullptr;
    vDSP_DFT_SetupD complexInverse = nullptr;
    FFTSetupD fftSetup = nullptr;
    double planningSeconds = 0;
    std::vector<double> splitReal;
    std::vector<double> splitImag;
    std::vector<double> lineInputReal;
    std::vector<double> lineInputImag;
    std::vector<double> lineOutputReal;
    std::vector<double> lineOutputImag;
};

size_t product(const ArrayDimensions& dimensions) {
    return std::accumulate(dimensions.begin(),dimensions.end(),size_t{1},std::multiplies<size_t>());
}

std::vector<size_t> strides(const ArrayDimensions& dimensions) {
    std::vector<size_t> values(dimensions.size(),1);
    for (size_t i = 1; i < dimensions.size(); ++i) {
        values[i] = values[i-1]*dimensions[i-1];
    }
    return values;
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
    output[transformDimensions.back()] = realDimensions[transformDimensions.back()]/2 + 1;
    return output;
}

bool isPowerOfTwo(size_t value) {
    return value > 0 && (value & (value-1)) == 0;
}

vDSP_Length exactLog2(size_t value) {
    if (!isPowerOfTwo(value)) {
        throw std::invalid_argument("The vDSP legacy 2-D strategy requires power-of-two transform lengths.");
    }
    vDSP_Length exponent = 0;
    while ((size_t{1} << exponent) < value) {
        ++exponent;
    }
    return exponent;
}

size_t baseIndex(size_t ordinal, const ArrayDimensions& dimensions, const std::vector<size_t>& dimensionStrides, const std::vector<size_t>& excludedDimensions) {
    size_t base = 0;
    for (size_t dimension = 0; dimension < dimensions.size(); ++dimension) {
        if (std::find(excludedDimensions.begin(),excludedDimensions.end(),dimension) != excludedDimensions.end()) {
            continue;
        }
        const size_t coordinate = ordinal % dimensions[dimension];
        ordinal /= dimensions[dimension];
        base += coordinate*dimensionStrides[dimension];
    }
    return base;
}

size_t batchCount(const ArrayDimensions& dimensions, const std::vector<size_t>& excludedDimensions) {
    size_t count = 1;
    for (size_t dimension = 0; dimension < dimensions.size(); ++dimension) {
        if (std::find(excludedDimensions.begin(),excludedDimensions.end(),dimension) == excludedDimensions.end()) {
            count *= dimensions[dimension];
        }
    }
    return count;
}

double elapsedSeconds(const Clock::time_point& start, const Clock::time_point& end) {
    return std::chrono::duration<double>(end-start).count();
}

void destroyPlan(PlanHandle* handle) {
    if (!handle) return;
    if (handle->realForward) vDSP_DFT_DestroySetupD(handle->realForward);
    if (handle->realInverse) vDSP_DFT_DestroySetupD(handle->realInverse);
    if (handle->complexForward) vDSP_DFT_DestroySetupD(handle->complexForward);
    if (handle->complexInverse) vDSP_DFT_DestroySetupD(handle->complexInverse);
    if (handle->fftSetup) vDSP_destroy_fftsetupD(handle->fftSetup);
    delete handle;
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
        if (!condition) fail(identifier,message);
    }

    PlanHandle* handleFrom(const Array& input) {
        const uint64_t value = static_cast<uint64_t>(input[0]);
        require(value != 0,"VDSPEngineBenchmark:InvalidPlan","The vDSP benchmark plan handle is invalid.");
        return reinterpret_cast<PlanHandle*>(value);
    }

    void validateDimensions(const ArrayDimensions& actual, const ArrayDimensions& expected, const std::string& label) {
        require(actual == expected,"VDSPEngineBenchmark:DimensionMismatch",label + " dimensions do not match the plan.");
    }

    void configureWorkBuffers(PlanHandle& handle) {
        const size_t nComplex = product(handle.complexDimensions);
        const size_t realLength = handle.realDimensions[handle.transformDimensions.back()];
        const size_t complexLength = handle.realDimensions[handle.transformDimensions.front()];
        const size_t maximumLine = std::max(realLength,complexLength);
        handle.splitReal.resize(nComplex);
        handle.splitImag.resize(nComplex);
        handle.lineInputReal.resize(maximumLine);
        handle.lineInputImag.resize(maximumLine);
        handle.lineOutputReal.resize(maximumLine);
        handle.lineOutputImag.resize(maximumLine);
    }

    void create(ArgumentList outputs, ArgumentList inputs) {
        require(inputs.size() == 4,"VDSPEngineBenchmark:InvalidInputCount","create expects dimensions, ordered transform dimensions, and strategy.");
        require(outputs.size() >= 1 && outputs.size() <= 5,"VDSPEngineBenchmark:InvalidOutputCount","create returns up to five outputs.");
        auto handle = std::make_unique<PlanHandle>();
        const auto dimensionValues = numericVector(inputs[1]);
        handle->realDimensions.assign(dimensionValues.begin(),dimensionValues.end());
        handle->transformDimensions = numericVector(inputs[2]);
        require(handle->transformDimensions.size() == 2 && handle->transformDimensions[0] != handle->transformDimensions[1],"VDSPEngineBenchmark:InvalidTransformDimensions","Exactly two distinct ordered transform dimensions are required.");
        for (size_t& dimension : handle->transformDimensions) {
            require(dimension <= handle->realDimensions.size(),"VDSPEngineBenchmark:InvalidTransformDimension","A transform dimension exceeds the input rank.");
            --dimension;
        }
        handle->complexDimensions = outputDimensions(handle->realDimensions,handle->transformDimensions);
        handle->strategy = textValue(inputs[3]);
        require(handle->strategy == "vdsp-staged" || handle->strategy == "vdsp-2d","VDSPEngineBenchmark:UnknownStrategy","Strategy must be vdsp-staged or vdsp-2d.");
        const size_t realLength = handle->realDimensions[handle->transformDimensions.back()];
        const size_t complexLength = handle->realDimensions[handle->transformDimensions.front()];
        require(realLength % 2 == 0,"VDSPEngineBenchmark:OddRealLength","The compressed vDSP transform length must be even.");

        const auto planningStart = Clock::now();
        if (handle->strategy == "vdsp-staged") {
            handle->realForward = vDSP_DFT_zrop_CreateSetupD(nullptr,realLength,vDSP_DFT_FORWARD);
            handle->realInverse = vDSP_DFT_zrop_CreateSetupD(handle->realForward,realLength,vDSP_DFT_INVERSE);
            handle->complexForward = vDSP_DFT_zop_CreateSetupD(handle->realForward,complexLength,vDSP_DFT_FORWARD);
            handle->complexInverse = vDSP_DFT_zop_CreateSetupD(handle->realForward,complexLength,vDSP_DFT_INVERSE);
            require(handle->realForward && handle->realInverse && handle->complexForward && handle->complexInverse,"VDSPEngineBenchmark:SetupCreationFailed","Accelerate does not provide a DFT setup for this workload.");
        } else {
            const vDSP_Length realLog2 = exactLog2(realLength);
            const vDSP_Length complexLog2 = exactLog2(complexLength);
            handle->fftSetup = vDSP_create_fftsetupD(std::max(realLog2,complexLog2),kFFTRadix2);
            require(handle->fftSetup != nullptr,"VDSPEngineBenchmark:SetupCreationFailed","Accelerate failed to create the legacy 2-D FFT setup.");
        }
        const auto planningEnd = Clock::now();
        handle->planningSeconds = elapsedSeconds(planningStart,planningEnd);
        configureWorkBuffers(*handle);

        if (outputs.size() > 0) outputs[0] = factory.createScalar(reinterpret_cast<uint64_t>(handle.release()));
        PlanHandle* returnedHandle = reinterpret_cast<PlanHandle*>(static_cast<uint64_t>(outputs[0][0]));
        if (outputs.size() > 1) {
            auto dimensions = factory.createArray<double>({1,returnedHandle->complexDimensions.size()});
            std::copy(returnedHandle->complexDimensions.begin(),returnedHandle->complexDimensions.end(),dimensions.begin());
            outputs[1] = std::move(dimensions);
        }
        if (outputs.size() > 2) outputs[2] = factory.createScalar(1.0/static_cast<double>(realLength*complexLength));
        if (outputs.size() > 3) outputs[3] = factory.createScalar(returnedHandle->planningSeconds);
        if (outputs.size() > 4) outputs[4] = factory.createScalar(std::string("engine-managed"));
    }

    void free(ArgumentList inputs) {
        require(inputs.size() == 2,"VDSPEngineBenchmark:InvalidInputCount","free expects one plan handle.");
        destroyPlan(handleFrom(inputs[1]));
    }

    void stagedForward(PlanHandle& handle, const double* input, std::complex<double>* output, double& kernelSeconds) {
        const size_t realDimension = handle.transformDimensions.back();
        const size_t complexDimension = handle.transformDimensions.front();
        const size_t realLength = handle.realDimensions[realDimension];
        const size_t halfLength = realLength/2;
        const size_t complexLength = handle.realDimensions[complexDimension];
        const auto realStrides = strides(handle.realDimensions);
        const auto complexStrides = strides(handle.complexDimensions);

        const size_t realBatches = batchCount(handle.realDimensions,{realDimension});
        for (size_t batch = 0; batch < realBatches; ++batch) {
            const size_t realBase = baseIndex(batch,handle.realDimensions,realStrides,{realDimension});
            const size_t complexBase = baseIndex(batch,handle.complexDimensions,complexStrides,{realDimension});
            for (size_t index = 0; index < halfLength; ++index) {
                handle.lineInputReal[index] = input[realBase + (2*index)*realStrides[realDimension]];
                handle.lineInputImag[index] = input[realBase + (2*index+1)*realStrides[realDimension]];
            }
            const auto start = Clock::now();
            vDSP_DFT_ExecuteD(handle.realForward,handle.lineInputReal.data(),handle.lineInputImag.data(),handle.lineOutputReal.data(),handle.lineOutputImag.data());
            const auto end = Clock::now();
            kernelSeconds += elapsedSeconds(start,end);
            handle.splitReal[complexBase] = 0.5*handle.lineOutputReal[0];
            handle.splitImag[complexBase] = 0;
            for (size_t index = 1; index < halfLength; ++index) {
                const size_t outputIndex = complexBase + index*complexStrides[realDimension];
                handle.splitReal[outputIndex] = 0.5*handle.lineOutputReal[index];
                handle.splitImag[outputIndex] = 0.5*handle.lineOutputImag[index];
            }
            const size_t nyquistIndex = complexBase + halfLength*complexStrides[realDimension];
            handle.splitReal[nyquistIndex] = 0.5*handle.lineOutputImag[0];
            handle.splitImag[nyquistIndex] = 0;
        }

        const size_t complexBatches = batchCount(handle.complexDimensions,{complexDimension});
        for (size_t batch = 0; batch < complexBatches; ++batch) {
            const size_t base = baseIndex(batch,handle.complexDimensions,complexStrides,{complexDimension});
            for (size_t index = 0; index < complexLength; ++index) {
                const size_t splitIndex = base + index*complexStrides[complexDimension];
                handle.lineInputReal[index] = handle.splitReal[splitIndex];
                handle.lineInputImag[index] = handle.splitImag[splitIndex];
            }
            const auto start = Clock::now();
            vDSP_DFT_ExecuteD(handle.complexForward,handle.lineInputReal.data(),handle.lineInputImag.data(),handle.lineOutputReal.data(),handle.lineOutputImag.data());
            const auto end = Clock::now();
            kernelSeconds += elapsedSeconds(start,end);
            for (size_t index = 0; index < complexLength; ++index) {
                const size_t splitIndex = base + index*complexStrides[complexDimension];
                handle.splitReal[splitIndex] = handle.lineOutputReal[index];
                handle.splitImag[splitIndex] = handle.lineOutputImag[index];
            }
        }
        for (size_t index = 0; index < handle.splitReal.size(); ++index) {
            output[index] = std::complex<double>(handle.splitReal[index],handle.splitImag[index]);
        }
    }

    void stagedInverse(PlanHandle& handle, std::complex<double>* spectrum, double* output, double& kernelSeconds) {
        const size_t realDimension = handle.transformDimensions.back();
        const size_t complexDimension = handle.transformDimensions.front();
        const size_t realLength = handle.realDimensions[realDimension];
        const size_t halfLength = realLength/2;
        const size_t complexLength = handle.realDimensions[complexDimension];
        const auto realStrides = strides(handle.realDimensions);
        const auto complexStrides = strides(handle.complexDimensions);
        for (size_t index = 0; index < handle.splitReal.size(); ++index) {
            handle.splitReal[index] = spectrum[index].real();
            handle.splitImag[index] = spectrum[index].imag();
        }

        const size_t complexBatches = batchCount(handle.complexDimensions,{complexDimension});
        for (size_t batch = 0; batch < complexBatches; ++batch) {
            const size_t base = baseIndex(batch,handle.complexDimensions,complexStrides,{complexDimension});
            for (size_t index = 0; index < complexLength; ++index) {
                const size_t splitIndex = base + index*complexStrides[complexDimension];
                handle.lineInputReal[index] = handle.splitReal[splitIndex];
                handle.lineInputImag[index] = handle.splitImag[splitIndex];
            }
            const auto start = Clock::now();
            vDSP_DFT_ExecuteD(handle.complexInverse,handle.lineInputReal.data(),handle.lineInputImag.data(),handle.lineOutputReal.data(),handle.lineOutputImag.data());
            const auto end = Clock::now();
            kernelSeconds += elapsedSeconds(start,end);
            for (size_t index = 0; index < complexLength; ++index) {
                const size_t splitIndex = base + index*complexStrides[complexDimension];
                handle.splitReal[splitIndex] = handle.lineOutputReal[index];
                handle.splitImag[splitIndex] = handle.lineOutputImag[index];
            }
        }

        const size_t realBatches = batchCount(handle.realDimensions,{realDimension});
        for (size_t batch = 0; batch < realBatches; ++batch) {
            const size_t realBase = baseIndex(batch,handle.realDimensions,realStrides,{realDimension});
            const size_t complexBase = baseIndex(batch,handle.complexDimensions,complexStrides,{realDimension});
            handle.lineInputReal[0] = handle.splitReal[complexBase];
            const size_t nyquistIndex = complexBase + halfLength*complexStrides[realDimension];
            handle.lineInputImag[0] = handle.splitReal[nyquistIndex];
            for (size_t index = 1; index < halfLength; ++index) {
                const size_t splitIndex = complexBase + index*complexStrides[realDimension];
                handle.lineInputReal[index] = handle.splitReal[splitIndex];
                handle.lineInputImag[index] = handle.splitImag[splitIndex];
            }
            const auto start = Clock::now();
            vDSP_DFT_ExecuteD(handle.realInverse,handle.lineInputReal.data(),handle.lineInputImag.data(),handle.lineOutputReal.data(),handle.lineOutputImag.data());
            const auto end = Clock::now();
            kernelSeconds += elapsedSeconds(start,end);
            for (size_t index = 0; index < halfLength; ++index) {
                output[realBase + (2*index)*realStrides[realDimension]] = handle.lineOutputReal[index];
                output[realBase + (2*index+1)*realStrides[realDimension]] = handle.lineOutputImag[index];
            }
        }
        spectrum[0] = std::complex<double>(std::numeric_limits<double>::quiet_NaN(),0);
    }

    void twoDimensionalForward(PlanHandle& handle, const double* input, std::complex<double>* output, double& kernelSeconds) {
        const size_t realDimension = handle.transformDimensions.back();
        const size_t complexDimension = handle.transformDimensions.front();
        const size_t realLength = handle.realDimensions[realDimension];
        const size_t halfLength = realLength/2;
        const size_t complexLength = handle.realDimensions[complexDimension];
        const auto realStrides = strides(handle.realDimensions);
        const auto complexStrides = strides(handle.complexDimensions);
        const size_t planeSize = halfLength*complexLength;
        const vDSP_Length realLog2 = exactLog2(realLength);
        const vDSP_Length complexLog2 = exactLog2(complexLength);
        DSPDoubleSplitComplex plane{handle.splitReal.data(),handle.splitImag.data()};
        const size_t batches = batchCount(handle.realDimensions,{realDimension,complexDimension});
        for (size_t batch = 0; batch < batches; ++batch) {
            const size_t realBase = baseIndex(batch,handle.realDimensions,realStrides,{realDimension,complexDimension});
            const size_t complexBase = baseIndex(batch,handle.complexDimensions,complexStrides,{realDimension,complexDimension});
            for (size_t complexIndex = 0; complexIndex < complexLength; ++complexIndex) {
                for (size_t realIndex = 0; realIndex < halfLength; ++realIndex) {
                    const size_t planeIndex = complexIndex*halfLength + realIndex;
                    handle.splitReal[planeIndex] = input[realBase + complexIndex*realStrides[complexDimension] + (2*realIndex)*realStrides[realDimension]];
                    handle.splitImag[planeIndex] = input[realBase + complexIndex*realStrides[complexDimension] + (2*realIndex+1)*realStrides[realDimension]];
                }
            }
            const auto start = Clock::now();
            vDSP_fft2d_zripD(handle.fftSetup,&plane,1,static_cast<vDSP_Stride>(halfLength),realLog2,complexLog2,FFT_FORWARD);
            const auto end = Clock::now();
            kernelSeconds += elapsedSeconds(start,end);

            for (size_t complexIndex = 0; complexIndex < complexLength; ++complexIndex) {
                for (size_t realIndex = 1; realIndex < halfLength; ++realIndex) {
                    const size_t planeIndex = complexIndex*halfLength + realIndex;
                    const size_t outputIndex = complexBase + complexIndex*complexStrides[complexDimension] + realIndex*complexStrides[realDimension];
                    output[outputIndex] = 0.5*std::complex<double>(handle.splitReal[planeIndex],handle.splitImag[planeIndex]);
                }
            }
            output[complexBase] = std::complex<double>(0.5*handle.splitReal[0],0);
            output[complexBase + halfLength*complexStrides[realDimension]] = std::complex<double>(0.5*handle.splitImag[0],0);
            const size_t complexNyquist = complexLength/2;
            output[complexBase + complexNyquist*complexStrides[complexDimension]] = std::complex<double>(0.5*handle.splitReal[halfLength],0);
            output[complexBase + complexNyquist*complexStrides[complexDimension] + halfLength*complexStrides[realDimension]] = std::complex<double>(0.5*handle.splitImag[halfLength],0);
            for (size_t complexIndex = 1; complexIndex < complexNyquist; ++complexIndex) {
                const size_t firstPacked = (2*complexIndex)*halfLength;
                const size_t secondPacked = (2*complexIndex+1)*halfLength;
                const std::complex<double> zeroValue = 0.5*std::complex<double>(handle.splitReal[firstPacked],handle.splitReal[secondPacked]);
                const std::complex<double> nyquistValue = 0.5*std::complex<double>(handle.splitImag[firstPacked],handle.splitImag[secondPacked]);
                output[complexBase + complexIndex*complexStrides[complexDimension]] = zeroValue;
                output[complexBase + (complexLength-complexIndex)*complexStrides[complexDimension]] = std::conj(zeroValue);
                output[complexBase + complexIndex*complexStrides[complexDimension] + halfLength*complexStrides[realDimension]] = nyquistValue;
                output[complexBase + (complexLength-complexIndex)*complexStrides[complexDimension] + halfLength*complexStrides[realDimension]] = std::conj(nyquistValue);
            }
            std::fill(handle.splitReal.begin(),handle.splitReal.begin()+planeSize,0);
            std::fill(handle.splitImag.begin(),handle.splitImag.begin()+planeSize,0);
        }
    }

    void twoDimensionalInverse(PlanHandle& handle, std::complex<double>* spectrum, double* output, double& kernelSeconds) {
        const size_t realDimension = handle.transformDimensions.back();
        const size_t complexDimension = handle.transformDimensions.front();
        const size_t realLength = handle.realDimensions[realDimension];
        const size_t halfLength = realLength/2;
        const size_t complexLength = handle.realDimensions[complexDimension];
        const auto realStrides = strides(handle.realDimensions);
        const auto complexStrides = strides(handle.complexDimensions);
        const vDSP_Length realLog2 = exactLog2(realLength);
        const vDSP_Length complexLog2 = exactLog2(complexLength);
        DSPDoubleSplitComplex plane{handle.splitReal.data(),handle.splitImag.data()};
        const size_t batches = batchCount(handle.realDimensions,{realDimension,complexDimension});
        for (size_t batch = 0; batch < batches; ++batch) {
            const size_t realBase = baseIndex(batch,handle.realDimensions,realStrides,{realDimension,complexDimension});
            const size_t complexBase = baseIndex(batch,handle.complexDimensions,complexStrides,{realDimension,complexDimension});
            for (size_t complexIndex = 0; complexIndex < complexLength; ++complexIndex) {
                for (size_t realIndex = 1; realIndex < halfLength; ++realIndex) {
                    const size_t planeIndex = complexIndex*halfLength + realIndex;
                    const std::complex<double> value = spectrum[complexBase + complexIndex*complexStrides[complexDimension] + realIndex*complexStrides[realDimension]];
                    handle.splitReal[planeIndex] = value.real();
                    handle.splitImag[planeIndex] = value.imag();
                }
            }
            handle.splitReal[0] = spectrum[complexBase].real();
            handle.splitImag[0] = spectrum[complexBase + halfLength*complexStrides[realDimension]].real();
            const size_t complexNyquist = complexLength/2;
            handle.splitReal[halfLength] = spectrum[complexBase + complexNyquist*complexStrides[complexDimension]].real();
            handle.splitImag[halfLength] = spectrum[complexBase + complexNyquist*complexStrides[complexDimension] + halfLength*complexStrides[realDimension]].real();
            for (size_t complexIndex = 1; complexIndex < complexNyquist; ++complexIndex) {
                const std::complex<double> zeroValue = spectrum[complexBase + complexIndex*complexStrides[complexDimension]];
                const std::complex<double> nyquistValue = spectrum[complexBase + complexIndex*complexStrides[complexDimension] + halfLength*complexStrides[realDimension]];
                const size_t firstPacked = (2*complexIndex)*halfLength;
                const size_t secondPacked = (2*complexIndex+1)*halfLength;
                handle.splitReal[firstPacked] = zeroValue.real();
                handle.splitReal[secondPacked] = zeroValue.imag();
                handle.splitImag[firstPacked] = nyquistValue.real();
                handle.splitImag[secondPacked] = nyquistValue.imag();
            }
            const auto start = Clock::now();
            vDSP_fft2d_zripD(handle.fftSetup,&plane,1,static_cast<vDSP_Stride>(halfLength),realLog2,complexLog2,FFT_INVERSE);
            const auto end = Clock::now();
            kernelSeconds += elapsedSeconds(start,end);
            for (size_t complexIndex = 0; complexIndex < complexLength; ++complexIndex) {
                for (size_t realIndex = 0; realIndex < halfLength; ++realIndex) {
                    const size_t planeIndex = complexIndex*halfLength + realIndex;
                    output[realBase + complexIndex*realStrides[complexDimension] + (2*realIndex)*realStrides[realDimension]] = handle.splitReal[planeIndex];
                    output[realBase + complexIndex*realStrides[complexDimension] + (2*realIndex+1)*realStrides[realDimension]] = handle.splitImag[planeIndex];
                }
            }
        }
        spectrum[0] = std::complex<double>(std::numeric_limits<double>::quiet_NaN(),0);
    }

    void r2c(ArgumentList outputs, ArgumentList inputs) {
        require(inputs.size() == 3 || inputs.size() == 4,"VDSPEngineBenchmark:InvalidInputCount","r2c expects a plan, real input, and optional preallocated spectrum.");
        require(outputs.size() >= 1 && outputs.size() <= 5,"VDSPEngineBenchmark:InvalidOutputCount","r2c returns the spectrum and timing/alignment metadata.");
        PlanHandle* handle = handleFrom(inputs[1]);
        require(inputs[2].getType() == matlab::data::ArrayType::DOUBLE,"VDSPEngineBenchmark:InvalidRealInput","r2c input must contain real doubles.");
        validateDimensions(inputs[2].getDimensions(),handle->realDimensions,"Real input");
        const double* input = readPointer<double>(inputs[2]);
        TypedArray<std::complex<double>> output = inputs.size() == 4 ? TypedArray<std::complex<double>>(std::move(inputs[3])) : factory.createArray<std::complex<double>>(handle->complexDimensions);
        validateDimensions(output.getDimensions(),handle->complexDimensions,"Spectrum output");
        std::complex<double>* outputPointer = &(*output.begin());
        double kernelSeconds = 0;
        const auto pipelineStart = Clock::now();
        if (handle->strategy == "vdsp-staged") stagedForward(*handle,input,outputPointer,kernelSeconds);
        else twoDimensionalForward(*handle,input,outputPointer,kernelSeconds);
        const auto pipelineEnd = Clock::now();
        outputs[0] = std::move(output);
        if (outputs.size() > 1) outputs[1] = factory.createScalar(kernelSeconds);
        if (outputs.size() > 2) outputs[2] = factory.createScalar(elapsedSeconds(pipelineStart,pipelineEnd));
        if (outputs.size() > 3) outputs[3] = factory.createScalar(static_cast<double>(reinterpret_cast<uintptr_t>(input) % 16));
        if (outputs.size() > 4) outputs[4] = factory.createScalar(static_cast<double>(reinterpret_cast<uintptr_t>(outputPointer) % 16));
    }

    void c2r(ArgumentList outputs, ArgumentList inputs) {
        require(inputs.size() == 4,"VDSPEngineBenchmark:InvalidInputCount","c2r expects a plan, destructive spectrum input, and preallocated real output.");
        require(outputs.size() >= 2 && outputs.size() <= 6,"VDSPEngineBenchmark:InvalidOutputCount","c2r returns the destroyed spectrum, real output, and timing/alignment metadata.");
        PlanHandle* handle = handleFrom(inputs[1]);
        TypedArray<std::complex<double>> spectrum = std::move(inputs[2]);
        TypedArray<double> output = std::move(inputs[3]);
        validateDimensions(spectrum.getDimensions(),handle->complexDimensions,"Spectrum input");
        validateDimensions(output.getDimensions(),handle->realDimensions,"Real output");
        std::complex<double>* spectrumPointer = &(*spectrum.begin());
        double* outputPointer = &(*output.begin());
        double kernelSeconds = 0;
        const auto pipelineStart = Clock::now();
        if (handle->strategy == "vdsp-staged") stagedInverse(*handle,spectrumPointer,outputPointer,kernelSeconds);
        else twoDimensionalInverse(*handle,spectrumPointer,outputPointer,kernelSeconds);
        const auto pipelineEnd = Clock::now();
        outputs[0] = std::move(spectrum);
        outputs[1] = std::move(output);
        if (outputs.size() > 2) outputs[2] = factory.createScalar(kernelSeconds);
        if (outputs.size() > 3) outputs[3] = factory.createScalar(elapsedSeconds(pipelineStart,pipelineEnd));
        if (outputs.size() > 4) outputs[4] = factory.createScalar(static_cast<double>(reinterpret_cast<uintptr_t>(spectrumPointer) % 16));
        if (outputs.size() > 5) outputs[5] = factory.createScalar(static_cast<double>(reinterpret_cast<uintptr_t>(outputPointer) % 16));
    }

    void info(ArgumentList outputs) {
        require(outputs.size() == 2,"VDSPEngineBenchmark:InvalidOutputCount","info returns the engine name and loaded framework path.");
        Dl_info information{};
        std::string libraryPath = "Accelerate.framework";
        if (dladdr(reinterpret_cast<void*>(vDSP_DFT_ExecuteD),&information) != 0 && information.dli_fname) libraryPath = information.dli_fname;
        outputs[0] = factory.createScalar(std::string("Apple Accelerate/vDSP"));
        outputs[1] = factory.createScalar(libraryPath);
    }

public:
    void operator()(ArgumentList outputs, ArgumentList inputs) override {
        require(!inputs.empty() && inputs[0].getType() == matlab::data::ArrayType::CHAR,"VDSPEngineBenchmark:InvalidCommand","The first input must be a command character vector.");
        const std::string command = textValue(inputs[0]);
        try {
            if (command == "create") create(outputs,inputs);
            else if (command == "free") free(inputs);
            else if (command == "r2c") r2c(outputs,inputs);
            else if (command == "c2r") c2r(outputs,inputs);
            else if (command == "info") info(outputs);
            else fail("VDSPEngineBenchmark:UnknownCommand","Unknown vDSP benchmark command: " + command);
        } catch (const matlab::engine::MATLABException&) {
            throw;
        } catch (const std::exception& exception) {
            fail("VDSPEngineBenchmark:NativeFailure",exception.what());
        }
    }
};
