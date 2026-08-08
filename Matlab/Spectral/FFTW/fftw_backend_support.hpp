#pragma once

#include "mex.hpp"
#include <fftw3.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

namespace fftw_backend {

using Clock = std::chrono::steady_clock;
using matlab::data::Array;
using matlab::data::ArrayDimensions;
using matlab::data::CharArray;
using matlab::data::TypedArray;

struct ScratchBuffer {
    double* base = nullptr;
    double* data = nullptr;
};

inline size_t product(const ArrayDimensions& dimensions) {
    return std::accumulate(dimensions.begin(),dimensions.end(),size_t{1},std::multiplies<size_t>());
}

inline ArrayDimensions canonicalDimensions(ArrayDimensions dimensions) {
    while (dimensions.size() > 2 && dimensions.back() == 1) dimensions.pop_back();
    return dimensions;
}

inline bool dimensionsEqualIgnoringTrailingSingletons(const ArrayDimensions& first, const ArrayDimensions& second) {
    return canonicalDimensions(first) == canonicalDimensions(second);
}

inline std::vector<size_t> strides(const ArrayDimensions& dimensions) {
    std::vector<size_t> values(dimensions.size(),1);
    for (size_t i = 1; i < dimensions.size(); ++i) values[i] = values[i-1]*dimensions[i-1];
    return values;
}

inline int checkedInt(size_t value, const char* label) {
    if (value > static_cast<size_t>(std::numeric_limits<int>::max())) throw std::invalid_argument(std::string(label) + " exceeds the FFTW guru integer range.");
    return static_cast<int>(value);
}

template <typename T>
const T* readPointer(const Array& input) {
    const TypedArray<T> typed = input;
    matlab::data::TypedIterator<const T> iterator(typed.begin());
    return iterator.operator->();
}

inline std::vector<size_t> numericVector(const Array& input) {
    TypedArray<double> values = input;
    std::vector<size_t> output;
    output.reserve(values.getNumberOfElements());
    for (double value : values) {
        if (!std::isfinite(value) || value < 1 || value != std::floor(value)) throw std::invalid_argument("Dimension values must be positive integers.");
        output.push_back(static_cast<size_t>(value));
    }
    return output;
}

inline std::string textValue(const Array& input) {
    CharArray value = input;
    return value.toAscii();
}

inline ArrayDimensions outputDimensions(const ArrayDimensions& realDimensions, const std::vector<size_t>& transformDimensions) {
    if (transformDimensions.empty()) throw std::invalid_argument("At least one transform dimension is required.");
    ArrayDimensions output = realDimensions;
    const size_t compressedDimension = transformDimensions.back();
    output[compressedDimension] = realDimensions[compressedDimension]/2 + 1;
    return output;
}

inline std::vector<fftw_iodim> transformIODimensions(const ArrayDimensions& realDimensions, const ArrayDimensions& complexDimensions, const std::vector<size_t>& transformDimensions, bool forward) {
    const auto realStrides = strides(realDimensions);
    const auto complexStrides = strides(complexDimensions);
    std::vector<fftw_iodim> dimensions;
    dimensions.reserve(transformDimensions.size());
    for (size_t dimension : transformDimensions) {
        dimensions.push_back(fftw_iodim{checkedInt(realDimensions[dimension],"Transform length"),checkedInt(forward ? realStrides[dimension] : complexStrides[dimension],"Input stride"),checkedInt(forward ? complexStrides[dimension] : realStrides[dimension],"Output stride")});
    }
    return dimensions;
}

inline std::vector<fftw_iodim> batchIODimensions(const ArrayDimensions& realDimensions, const ArrayDimensions& complexDimensions, const std::vector<size_t>& transformDimensions, bool forward) {
    const auto realStrides = strides(realDimensions);
    const auto complexStrides = strides(complexDimensions);
    std::vector<fftw_iodim> dimensions;
    for (size_t dimension = 0; dimension < realDimensions.size(); ++dimension) {
        if (std::find(transformDimensions.begin(),transformDimensions.end(),dimension) != transformDimensions.end()) continue;
        dimensions.push_back(fftw_iodim{checkedInt(realDimensions[dimension],"Batch length"),checkedInt(forward ? realStrides[dimension] : complexStrides[dimension],"Batch input stride"),checkedInt(forward ? complexStrides[dimension] : realStrides[dimension],"Batch output stride")});
    }
    return dimensions;
}

inline ScratchBuffer alignedBuffer(size_t nDoubles, int targetAlignment) {
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

inline ScratchBuffer alignedBufferWithOffset(size_t nDoubles, size_t offsetDoubles, int targetAlignment) {
    constexpr size_t padding = 64;
    ScratchBuffer buffer;
    buffer.base = fftw_alloc_real(nDoubles + padding);
    if (!buffer.base) throw std::bad_alloc();
    for (size_t offset = 0; offset < padding; ++offset) {
        if (fftw_alignment_of(buffer.base + offset + offsetDoubles) == targetAlignment) {
            buffer.data = buffer.base + offset;
            return buffer;
        }
    }
    fftw_free(buffer.base);
    throw std::runtime_error("Unable to manufacture the requested offset FFTW alignment class.");
}

inline double elapsedSeconds(const Clock::time_point& start, const Clock::time_point& end) {
    return std::chrono::duration<double>(end-start).count();
}

inline bool hasDistinctDimensions(const std::vector<size_t>& dimensions) {
    std::vector<size_t> sorted = dimensions;
    std::sort(sorted.begin(),sorted.end());
    return std::adjacent_find(sorted.begin(),sorted.end()) == sorted.end();
}

} // namespace fftw_backend
